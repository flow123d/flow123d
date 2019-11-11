/*
 * darcy_flow_assembly.hh
 *
 *  Created on: Apr 21, 2016
 *      Author: jb
 */

#ifndef SRC_FLOW_DARCY_FLOW_ASSEMBLY_LMH_HH_
#define SRC_FLOW_DARCY_FLOW_ASSEMBLY_LMH_HH_

#include "mesh/long_idx.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/neighbours.h"
#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_values_views.hh"
#include "fem/fe_system.hh"
#include "quadrature/quadrature_lib.hh"
#include "flow/mh_dofhandler.hh"

#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "la/linsys_BDDC.hh"
#include "la/schur.hh"

#include "la/local_system.hh"

#include "coupling/balance.hh"
#include "flow/assembly_mh.hh"
#include "flow/darcy_flow_lmh.hh"
#include "flow/mortar_assembly.hh"

/** Copy of the assembly class for MH implementation,
 * with Lumping and further improvements.
 * Used also for Richards.
 */
template<int dim>
class AssemblyLMH : public AssemblyBase
{
public:
    typedef std::shared_ptr<DarcyLMH::EqData> AssemblyDataPtrLMH;
    
    AssemblyLMH<dim>(AssemblyDataPtrLMH data)
    : quad_(dim, 3),
        fe_values_(map_, quad_, fe_rt_,
                update_values | update_gradients | update_JxW_values | update_quadrature_points),

        velocity_interpolation_quad_(dim, 0), // veloctiy values in barycenter
        velocity_interpolation_fv_(map_,velocity_interpolation_quad_, fe_rt_, update_values | update_quadrature_points),

        ad_(data),
        loc_system_(size(), size()),
        loc_system_vb_(2,2),
        edge_indices_(dim+1)
    {
        // local numbering of dofs for MH system
        // note: this shortcut supposes that the fe_system is the same on all elements
        // the function DiscreteSpace.fe(ElementAccessor) does not in fact depend on the element accessor
        auto fe = ad_->dh_->ds()->fe(ad_->dh_->own_range().begin()->elm()).get<dim>();
        FESystem<dim>* fe_system = dynamic_cast<FESystem<dim>*>(fe.get());
        loc_side_dofs = fe_system->fe_dofs(0);;
        loc_ele_dof = fe_system->fe_dofs(1)[0];
        loc_edge_dofs = fe_system->fe_dofs(2);;
        
        unsigned int nsides = dim+1;
        
        // create local sparsity pattern
        arma::umat sp(size(),size());
        sp.zeros();
        sp.submat(0, 0, nsides, nsides).ones();
        sp.diag().ones();
        // armadillo 8.4.3 bug with negative sub diagonal index
        // sp.diag(nsides+1).ones();
        // sp.diag(-nsides-1).ones();
        // sp.print();
        
        sp.submat(0, nsides+1, nsides-1, size()-1).diag().ones();
        sp.submat(nsides+1, 0, size()-1, nsides-1).diag().ones();
        
        loc_system_.set_sparsity(sp);
        
        // local system 2x2 for vb neighbourings is full matrix
        // this matrix cannot be influenced by any BC (no elimination can take place)
        sp.ones(2,2);
        loc_system_vb_.set_sparsity(sp);

        if (ad_->mortar_method_ == DarcyFlowInterface::MortarP0) {
            mortar_assembly = std::make_shared<P0_CouplingAssembler>(ad_);
        } else if (ad_->mortar_method_ == DarcyFlowInterface::MortarP1) {
            mortar_assembly = std::make_shared<P1_CouplingAssembler>(ad_);
        }

    }


    ~AssemblyLMH<dim>() override
    {}

//     LocalSystem& get_local_system() override
//         { return loc_system_;}
    
    void fix_velocity(LocalElementAccessorBase<3> ele_ac) override
    {
        if (mortar_assembly)
            mortar_assembly->fix_velocity(ele_ac);
    }

    void assemble(LocalElementAccessorBase<3> ele_ac) override
    {
        ASSERT_EQ_DBG(ele_ac.dim(), dim);
        loc_system_.reset();
    
        set_dofs_and_bc(ele_ac);
        
        assemble_sides(ele_ac);
        assemble_element(ele_ac);
        assemble_source_term(ele_ac);
        
        ad_->lin_sys->set_local_system(loc_system_);

        assembly_dim_connections(ele_ac);

        if (ad_->balance != nullptr)
            add_fluxes_in_balance_matrix(ele_ac);

        if (mortar_assembly)
            mortar_assembly->assembly(ele_ac);
    }

    void assembly_local_vb(ElementAccessor<3> ele, DHCellSide neighb_side) override
    {
        ASSERT_LT_DBG(ele->dim(), 3);
        //DebugOut() << "alv " << print_var(this);
        //START_TIMER("Assembly<dim>::assembly_local_vb");
        // compute normal vector to side
        arma::vec3 nv;
        ElementAccessor<3> ele_higher = ad_->mesh->element_accessor( neighb_side.element().idx() );
        ngh_values_.fe_side_values_.reinit(ele_higher, neighb_side.side_idx());
        nv = ngh_values_.fe_side_values_.normal_vector(0);

        double value = ad_->sigma.value( ele.centre(), ele) *
                        2*ad_->conductivity.value( ele.centre(), ele) *
                        arma::dot(ad_->anisotropy.value( ele.centre(), ele)*nv, nv) *
                        ad_->cross_section.value( neighb_side.centre(), ele_higher ) * // cross-section of higher dim. (2d)
                        ad_->cross_section.value( neighb_side.centre(), ele_higher ) /
                        ad_->cross_section.value( ele.centre(), ele ) *      // crossection of lower dim.
						neighb_side.measure();

        loc_system_vb_.add_value(0,0, -value);
        loc_system_vb_.add_value(0,1,  value);
        loc_system_vb_.add_value(1,0,  value);
        loc_system_vb_.add_value(1,1, -value);
    }


    arma::vec3 make_element_vector(LocalElementAccessorBase<3> ele_ac) override
    {
        //START_TIMER("Assembly<dim>::make_element_vector");
        arma::vec3 flux_in_center;
        flux_in_center.zeros();
        auto ele = ele_ac.element_accessor();

        velocity_interpolation_fv_.reinit(ele);
        for (unsigned int li = 0; li < ele->n_sides(); li++) {
            flux_in_center += ad_->data_vec_[ ele_ac.side_local_row(li) ]
                        * velocity_interpolation_fv_.vector_view(0).value(li,0);
        }


        flux_in_center /= ad_->cross_section.value(ele.centre(), ele );
        return flux_in_center;
    }

    void postprocess_velocity(const DHCellAccessor& dh_cell)
    {
        ElementAccessor<3> ele = dh_cell.elm();
        
        double edge_scale = ele.measure()
                              * ad_->cross_section.value(ele.centre(), ele)
                              / ele->n_sides();
        
        double edge_source_term = edge_scale * ad_->water_source_density.value(ele.centre(), ele);
      
        postprocess_velocity_specific(dh_cell, edge_scale, edge_source_term);
    }
    
protected:
    
    static const unsigned int size()
    {
        // dofs: velocity, pressure, edge pressure
        return RefElement<dim>::n_sides + 1 + RefElement<dim>::n_sides;
    }
    
    void set_dofs_and_bc(LocalElementAccessorBase<3> ele_ac){
        
        ASSERT_DBG(ele_ac.dim() == dim);
        
        //set global dof for element (pressure)
        loc_system_.row_dofs[loc_ele_dof] = loc_system_.col_dofs[loc_ele_dof] = ele_ac.ele_row();
        
        //shortcuts
        const unsigned int nsides = ele_ac.n_sides();
        LinSys *ls = ad_->lin_sys;
        
        Boundary *bcd;
        unsigned int side_row, edge_row;
        
        dirichlet_edge.resize(nsides);
        for (unsigned int i = 0; i < nsides; i++) {

            side_row = loc_side_dofs[i];    //local
            edge_row = loc_edge_dofs[i];    //local
            loc_system_.row_dofs[side_row] = loc_system_.col_dofs[side_row] = ele_ac.side_row(i);    //global
            loc_system_.row_dofs[edge_row] = loc_system_.col_dofs[edge_row] = ele_ac.edge_row(i);    //global
            
            bcd = ele_ac.side(i)->cond();
            dirichlet_edge[i] = 0;
            if (bcd) {
                ElementAccessor<3> b_ele = bcd->element_accessor();
                DarcyMH::EqData::BC_Type type = (DarcyMH::EqData::BC_Type)ad_->bc_type.value(b_ele.centre(), b_ele);

                double cross_section = ad_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor());

                if ( type == DarcyMH::EqData::none) {
                    // homogeneous neumann
                } else if ( type == DarcyMH::EqData::dirichlet ) {
                    double bc_pressure = ad_->bc_pressure.value(b_ele.centre(), b_ele);
                    loc_system_.set_solution(loc_edge_dofs[i],bc_pressure,-1);
                    dirichlet_edge[i] = 1;
                    
                } else if ( type == DarcyMH::EqData::total_flux) {
                    // internally we work with outward flux
                    double bc_flux = -ad_->bc_flux.value(b_ele.centre(), b_ele);
                    double bc_pressure = ad_->bc_pressure.value(b_ele.centre(), b_ele);
                    double bc_sigma = ad_->bc_robin_sigma.value(b_ele.centre(), b_ele);
                    
                    dirichlet_edge[i] = 2;  // to be skipped in LMH source assembly
                    loc_system_.add_value(edge_row, edge_row,
                                            -b_ele.measure() * bc_sigma * cross_section,
                                            (bc_flux - bc_sigma * bc_pressure) * b_ele.measure() * cross_section);
                }
                else if (type==DarcyMH::EqData::seepage) {
                    ad_->is_linear=false;

                    unsigned int loc_edge_idx = bcd->bc_ele_idx_;
                    char & switch_dirichlet = ad_->bc_switch_dirichlet[loc_edge_idx];
                    double bc_pressure = ad_->bc_switch_pressure.value(b_ele.centre(), b_ele);
                    double bc_flux = -ad_->bc_flux.value(b_ele.centre(), b_ele);
                    double side_flux = bc_flux * b_ele.measure() * cross_section;

                    // ** Update BC type. **
                    if (switch_dirichlet) {
                        // check and possibly switch to flux BC
                        // The switch raise error on the corresponding edge row.
                        // Magnitude of the error is abs(solution_flux - side_flux).
                        ASSERT_DBG(ad_->dh_->distr()->is_local(ele_ac.side_row(i)))(ele_ac.side_row(i));
                        unsigned int loc_side_row = ele_ac.side_local_row(i);
                        double & solution_flux = ls->get_solution_array()[loc_side_row];

                        if ( solution_flux < side_flux) {
                            //DebugOut().fmt("x: {}, to neum, p: {} f: {} -> f: {}\n", b_ele.centre()[0], bc_pressure, solution_flux, side_flux);
                            solution_flux = side_flux;
                            switch_dirichlet=0;
                        }
                    } else {
                        // check and possibly switch to  pressure BC
                        // TODO: What is the appropriate DOF in not local?
                        // The switch raise error on the corresponding side row.
                        // Magnitude of the error is abs(solution_head - bc_pressure)
                        // Since usually K is very large, this error would be much
                        // higher then error caused by the inverse switch, this
                        // cause that a solution  with the flux violating the
                        // flux inequality leading may be accepted, while the error
                        // in pressure inequality is always satisfied.
                        ASSERT_DBG(ad_->dh_->distr()->is_local(ele_ac.edge_row(i)))(ele_ac.edge_row(i));
                        unsigned int loc_edge_row = ele_ac.edge_local_row(i);
                        double & solution_head = ls->get_solution_array()[loc_edge_row];

                        if ( solution_head > bc_pressure) {
                            //DebugOut().fmt("x: {}, to dirich, p: {} -> p: {} f: {}\n",b_ele.centre()[0], solution_head, bc_pressure, bc_flux);
                            solution_head = bc_pressure;
                            switch_dirichlet=1;
                        }
                    }
                    
                        // ** Apply BCUpdate BC type. **
                        // Force Dirichlet type during the first iteration of the unsteady case.
                        if (switch_dirichlet || ad_->force_bc_switch ) {
                            //DebugOut().fmt("x: {}, dirich, bcp: {}\n", b_ele.centre()[0], bc_pressure);
                            loc_system_.set_solution(loc_edge_dofs[i],bc_pressure, -1);
                            dirichlet_edge[i] = 1;
                        } else {
                            //DebugOut()("x: {}, neuman, q: {}  bcq: {}\n", b_ele.centre()[0], side_flux, bc_flux);
                            loc_system_.add_value(edge_row, side_flux);
                        }

                } else if (type==DarcyMH::EqData::river) {
                    ad_->is_linear=false;

                    double bc_pressure = ad_->bc_pressure.value(b_ele.centre(), b_ele);
                    double bc_switch_pressure = ad_->bc_switch_pressure.value(b_ele.centre(), b_ele);
                    double bc_flux = -ad_->bc_flux.value(b_ele.centre(), b_ele);
                    double bc_sigma = ad_->bc_robin_sigma.value(b_ele.centre(), b_ele);
                    ASSERT_DBG(ad_->dh_->distr()->is_local(ele_ac.edge_row(i)))(ele_ac.edge_row(i));
                    unsigned int loc_edge_row = ele_ac.edge_local_row(i);
                    double & solution_head = ls->get_solution_array()[loc_edge_row];

                    // Force Robin type during the first iteration of the unsteady case.
                    if (solution_head > bc_switch_pressure  || ad_->force_bc_switch) {
                        // Robin BC
                        //DebugOut().fmt("x: {}, robin, bcp: {}\n", b_ele.centre()[0], bc_pressure);
                        loc_system_.add_value(edge_row, edge_row,
                                                -b_ele.measure() * bc_sigma * cross_section,
												b_ele.measure() * cross_section * (bc_flux - bc_sigma * bc_pressure)  );
                    } else {
                        // Neumann BC
                        //DebugOut().fmt("x: {}, neuman, q: {}  bcq: {}\n", b_ele.centre()[0], bc_switch_pressure, bc_pressure);
                        double bc_total_flux = bc_flux + bc_sigma*(bc_switch_pressure - bc_pressure);
                        
                        loc_system_.add_value(edge_row, bc_total_flux * b_ele.measure() * cross_section);
                    }
                } 
                else {
                    xprintf(UsrErr, "BC type not supported.\n");
                }
            }
            loc_system_.add_value(side_row, edge_row, 1.0);
            loc_system_.add_value(edge_row, side_row, 1.0);
        }
        
//         DBGCOUT(<< "ele " << ele_ac.ele_global_idx() << ":  ");
//         for (unsigned int i = 0; i < nsides; i++) cout << loc_system_.row_dofs[loc_side_dofs[i]] << "  ";
//         cout << loc_system_.row_dofs[loc_ele_dof] << "  ";
//         for (unsigned int i = 0; i < nsides; i++) cout << loc_system_.row_dofs[loc_edge_dofs[i]] << "  ";
//         cout << "\n";
    }
        
     void assemble_sides(LocalElementAccessorBase<3> ele_ac) override
     {
        double cs = ad_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor());
        double conduct =  ad_->conductivity.value(ele_ac.centre(), ele_ac.element_accessor());
        double scale = 1 / cs /conduct;
        
        assemble_sides_scale(ele_ac, scale);
    }
    
    void assemble_sides_scale(LocalElementAccessorBase<3> ele_ac, double scale)
    {
        arma::vec3 &gravity_vec = ad_->gravity_vec_;
        
        ElementAccessor<3> ele =ele_ac.element_accessor();
        fe_values_.reinit(ele);
        unsigned int ndofs = fe_values_.get_fe()->n_dofs();
        unsigned int qsize = fe_values_.get_quadrature()->size();
        auto velocity = fe_values_.vector_view(0);

        for (unsigned int k=0; k<qsize; k++)
            for (unsigned int i=0; i<ndofs; i++){
                double rhs_val =
                        arma::dot(gravity_vec,velocity.value(i,k))
                        * fe_values_.JxW(k);
                loc_system_.add_value(i, rhs_val);
                
                for (unsigned int j=0; j<ndofs; j++){
                    double mat_val = 
                        arma::dot(velocity.value(i,k), //TODO: compute anisotropy before
                                    (ad_->anisotropy.value(ele_ac.centre(), ele_ac.element_accessor() )).i()
                                        * velocity.value(j,k))
                        * scale * fe_values_.JxW(k);
                    
                    loc_system_.add_value(i, j, mat_val);
                }
            }
        
        // assemble matrix for weights in BDDCML
        // approximation to diagonal of 
        // S = -C - B*inv(A)*B'
        // as 
        // diag(S) ~ - diag(C) - 1./diag(A)
        // the weights form a partition of unity to average a discontinuous solution from neighbouring subdomains
        // to a continuous one
        // it is important to scale the effect - if conductivity is low for one subdomain and high for the other,
        // trust more the one with low conductivity - it will be closer to the truth than an arithmetic average
        if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
            const arma::mat& local_matrix = loc_system_.get_matrix();
            for(unsigned int i=0; i < ndofs; i++) {
                double val_side =  local_matrix(i,i);
                double val_edge =  -1./local_matrix(i,i);

                unsigned int side_row = loc_system_.row_dofs[loc_side_dofs[i]];
                unsigned int edge_row = loc_system_.row_dofs[loc_edge_dofs[i]];
                static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( side_row, val_side );
                static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( edge_row, val_edge );
            }
        }
    }
    
    
    void assemble_element(LocalElementAccessorBase<3> ele_ac){
        // set block B, B': element-side, side-element
        
        for(unsigned int side = 0; side < loc_side_dofs.size(); side++){
            loc_system_.add_value(loc_ele_dof, loc_side_dofs[side], -1.0);
            loc_system_.add_value(loc_side_dofs[side], loc_ele_dof, -1.0);
        }
        
        if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
            double val_ele =  1.;
            static_cast<LinSys_BDDC*>(ad_->lin_sys)->
                            diagonal_weights_set_value( loc_system_.row_dofs[loc_ele_dof], val_ele );
        }
    }
    
    void assemble_source_term(LocalElementAccessorBase<3> ele_ac) override
    {
        ElementAccessor<3> ele =ele_ac.element_accessor();
        
        ele_ac.dh_cell().cell_with_other_dh(ad_->dh_cr_.get()).get_loc_dof_indices(edge_indices_);
        
        // compute lumped source
        double alpha = 1.0 / ele_ac.n_sides();
        double cross_section = ad_->cross_section.value(ele.centre(), ele);
        double coef = alpha * ele.measure() * cross_section;
        
        double source = ad_->water_source_density.value(ele_ac.centre(), ele_ac.element_accessor());
        double source_term = coef * source;
        
        // in unsteady, compute time term
        double storativity = 0.0;
        double time_term_diag = 0.0, time_term = 0.0, time_term_rhs = 0.0;
        
        if(! ad_->use_steady_assembly_)
        {
            storativity = ad_->storativity.value(ele_ac.centre(), ele_ac.element_accessor());
            time_term = coef * storativity;
        }
        
        for (unsigned int i=0; i<ele_ac.n_sides(); i++)
        {
            if(! ad_->use_steady_assembly_)
            {
                time_term_diag = time_term / ad_->time_step_;
                time_term_rhs = time_term_diag * ad_->previous_solution[ele_ac.edge_local_row(i)];
            }

            this->loc_system_.add_value(loc_edge_dofs[i], loc_edge_dofs[i],
                                        -time_term_diag,
                                        -source_term - time_term_rhs);
                
            if (ad_->balance != nullptr)
            {
                ad_->balance->add_source_values(ad_->water_balance_idx, ele.region().bulk_idx(), {(LongIdx)edge_indices_[i]}, {0},{source_term});
                if( ! ad_->use_steady_assembly_)
                {
                    ad_->balance->add_mass_matrix_values(ad_->water_balance_idx, ele_ac.region().bulk_idx(), { LongIdx(ele_ac.edge_row(i)) },
                                                         {time_term});
                }
            }
        }
    }

    void assembly_dim_connections(LocalElementAccessorBase<3> ele_ac){
        //D, E',E block: compatible connections: element-edge
        auto ele = ele_ac.element_accessor(); //ElementAccessor<3>
        DHCellAccessor dh_cell = ele_ac.dh_cell();
        
        // no Neighbours => nothing to asssemble here
        if(dh_cell.elm()->n_neighs_vb() == 0) return;
        
        int ele_row = ele_ac.ele_row();
        Neighbour *ngh;

        //DebugOut() << "adc " << print_var(this) << print_var(side_quad_.size());
        unsigned int i = 0;
        for ( DHCellSide neighb_side : dh_cell.neighb_sides() ) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            ngh = ele_ac.element_accessor()->neigh_vb[i];
            loc_system_vb_.reset();
            loc_system_vb_.row_dofs[0] = loc_system_vb_.col_dofs[0] = ele_row;
            DHCellAccessor cell_higher_dim = ele_ac.dh_cell().dh()->cell_accessor_from_element(neighb_side.elem_idx());
            LocalElementAccessorBase<3> acc_higher_dim( cell_higher_dim );
            for (unsigned int j = 0; j < neighb_side.element().dim()+1; j++)
            	if (neighb_side.element()->edge_idx(j) == ngh->edge_idx()) {
                    loc_system_vb_.row_dofs[1] = loc_system_vb_.col_dofs[1] = acc_higher_dim.edge_row(j);
            		break;
            	}

            assembly_local_vb(ele, neighb_side);

            ad_->lin_sys->set_local_system(loc_system_vb_);

            // update matrix for weights in BDDCML
            if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
               int ind = loc_system_vb_.row_dofs[1];
               // there is -value on diagonal in block C!
               double new_val = loc_system_vb_.get_matrix()(0,0);
               static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( ind, new_val );
            }
            ++i;
        }
    }

    void add_fluxes_in_balance_matrix(LocalElementAccessorBase<3> ele_ac){

        for (unsigned int i = 0; i < ele_ac.n_sides(); i++) {
            Boundary* bcd = ele_ac.side(i)->cond();

            if (bcd) {
                ad_->balance->add_flux_matrix_values(ad_->water_balance_idx, ele_ac.side(i),
                                                     {(LongIdx)(ele_ac.side_row(i))}, {1});
            }
        }
    }


    virtual void postprocess_velocity_specific(const DHCellAccessor& dh_cell, double edge_scale, double edge_source_term)// override
    {
        dh_cell.get_loc_dof_indices(indices_);
        ElementAccessor<3> ele = dh_cell.elm();
        
        double storativity = ad_->storativity.value(ele.centre(), ele);
        double new_pressure, old_pressure, time_term = 0.0;
        
        for (unsigned int i=0; i<ele->n_sides(); i++) {
            
            if( ! ad_->use_steady_assembly_)
            {
                new_pressure = ad_->data_vec_[         indices_[loc_edge_dofs[i]] ];
                old_pressure = ad_->previous_solution[ indices_[loc_edge_dofs[i]] ];
                time_term = edge_scale * storativity / ad_->time_step_ * (new_pressure - old_pressure);
            }
            
            ad_->data_vec_[indices_[loc_side_dofs[i]]] += edge_source_term - time_term;
        }
    }

    // assembly volume integrals
    FE_RT0<dim> fe_rt_;
    MappingP1<dim,3> map_;
    QGauss quad_;
    FEValues<dim,3> fe_values_;

    NeighSideValues<dim<3?dim:2> ngh_values_;

    // Interpolation of velocity into barycenters
    QGauss velocity_interpolation_quad_;
    FEValues<dim,3> velocity_interpolation_fv_;

    // data shared by assemblers of different dimension
    AssemblyDataPtrLMH ad_;
    
    /** TODO: Investigate why the hell do we need this flag.
    *  If removed, it does not break any of the integration tests,
    * however it must influence the Dirichlet rows in matrix.
    */
    std::vector<unsigned int> dirichlet_edge;

    LocalSystem loc_system_;
    LocalSystem loc_system_vb_;
    std::vector<unsigned int> loc_side_dofs;
    std::vector<unsigned int> loc_edge_dofs;
    unsigned int loc_ele_dof;

    std::shared_ptr<MortarAssemblyBase> mortar_assembly;
    
    // TODO: Update dofs only once, use the dofs from LocalSystem, once set_dofs and set_bc is separated.
    std::vector<int> indices_;
    std::vector<int> edge_indices_; ///< Dofs of discontinuous fields on element edges.
};


#endif /* SRC_FLOW_DARCY_FLOW_ASSEMBLY_LMH_HH_ */
