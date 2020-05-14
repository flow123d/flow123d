/*
 * darcy_flow_assembly.hh
 *
 *  Created on: Apr 21, 2016
 *      Author: jb
 */

#ifndef SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_
#define SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_

#include "system/index_types.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/neighbours.h"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_values_views.hh"
#include "quadrature/quadrature_lib.hh"

#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "la/linsys_BDDC.hh"
#include "la/schur.hh"

#include "la/local_system.hh"

#include "coupling/balance.hh"
#include "flow/darcy_flow_mh.hh"
#include "flow/mortar_assembly.hh"


/** Common abstract class for the assembly routines in Darcy flow. 
 * Is implemented in DarcyMH, DarcyLMH and RichardsLMH assembly classes,
 * which are independent of each other.
 */
class AssemblyBase
{
public:
    virtual void fix_velocity(const DHCellAccessor& dh_cell) = 0;
    virtual void assemble(const DHCellAccessor& dh_cell) = 0;
    virtual void assemble_reconstruct(const DHCellAccessor& dh_cell) = 0;

    /// Updates water content in Richards.
    virtual void update_water_content(const DHCellAccessor& dh_cell) = 0;

    /**
        * Generic creator of multidimensional assembly, i.e. vector of
        * particular assembly objects.
        */
    template< template<int dim> class Impl, class Data>
    static MultidimAssembly create(Data data) {
        return { std::make_shared<Impl<1> >(data),
            std::make_shared<Impl<2> >(data),
            std::make_shared<Impl<3> >(data) };
    }

    virtual ~AssemblyBase() {}
};


template <int dim>
class NeighSideValues {
private:
    // assembly face integrals (BC)
    QGauss side_quad_;
    FE_P_disc<dim+1> fe_p_disc_;
public:
    NeighSideValues<dim>()
    :  side_quad_(dim, 1),
       fe_p_disc_(0)
    {
        fe_side_values_.initialize(side_quad_, fe_p_disc_, update_normal_vectors);
    }
    FEValues<3> fe_side_values_;

};


/** MH version of Darcy flow assembly. It is supposed not to be improved anymore,
 * however it is kept functioning aside of the LMH lumped version until
 * the LMH version is stable and optimized.
 */
template<int dim>
class AssemblyMH : public AssemblyBase
{
public:
    typedef std::shared_ptr<DarcyMH::EqData>  AssemblyDataPtrMH;
    
    AssemblyMH<dim>(AssemblyDataPtrMH data)
    : quad_(dim, 3),
      velocity_interpolation_quad_(dim, 0), // veloctiy values in barycenter
      ad_(data),
      loc_system_(size(), size()),
      loc_system_vb_(2,2)

    {
        fe_values_.initialize(quad_, fe_rt_,
                update_values | update_gradients | update_JxW_values | update_quadrature_points);
        velocity_interpolation_fv_.initialize(velocity_interpolation_quad_, fe_rt_, update_values | update_quadrature_points);

        // local numbering of dofs for MH system
        unsigned int nsides = dim+1;
        loc_side_dofs.resize(nsides);
        loc_ele_dof = nsides;
        loc_edge_dofs.resize(nsides);
        for(unsigned int i = 0; i < nsides; i++){
            loc_side_dofs[i] = i;
            loc_edge_dofs[i] = nsides + i + 1;
        }
        //DebugOut() << print_var(this) << print_var(side_quad_.size());
        
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

    void assemble_reconstruct(const DHCellAccessor&) override
    {};
    void update_water_content(const DHCellAccessor&) override
    {};

    ~AssemblyMH<dim>() override
    {}

//     LocalSystem& get_local_system() override
//         { return loc_system_;}
    
    void fix_velocity(const DHCellAccessor& dh_cell) override
    {
        if (mortar_assembly)
            mortar_assembly->fix_velocity(dh_cell);
    }

    void assemble(const DHCellAccessor& dh_cell) override
    {
        ASSERT_EQ_DBG(dh_cell.dim(), dim);
        loc_system_.reset();
    
        set_dofs_and_bc(dh_cell);
        
        assemble_sides(dh_cell);
        assemble_element(dh_cell);
        
        loc_system_.eliminate_solution();
        ad_->lin_sys->set_local_system(loc_system_);

        assembly_dim_connections(dh_cell);

        if (ad_->balance != nullptr)
            add_fluxes_in_balance_matrix(dh_cell);

        if (mortar_assembly)
            mortar_assembly->assembly(dh_cell);
    }

    void assembly_local_vb(ElementAccessor<3> ele, DHCellSide neighb_side) //override
    {
        ASSERT_LT_DBG(ele->dim(), 3);
        //DebugOut() << "alv " << print_var(this);
        //START_TIMER("Assembly<dim>::assembly_local_vb");
        // compute normal vector to side
        arma::vec3 nv;
        ElementAccessor<3> ele_higher = ad_->mesh->element_accessor( neighb_side.element().idx() );
        ngh_values_.fe_side_values_.reinit(neighb_side.side());
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

protected:
    static unsigned int size()
    {
        // dofs: velocity, pressure, edge pressure
        return RefElement<dim>::n_sides + 1 + RefElement<dim>::n_sides;
    }
    
    void set_dofs_and_bc(const DHCellAccessor& dh_cell){

        local_dofs_ = dh_cell.get_loc_dof_indices();
        global_dofs_.resize(dh_cell.n_dofs());
        dh_cell.get_dof_indices(global_dofs_);

        const ElementAccessor<3> ele = dh_cell.elm();
        
        //set global dof for element (pressure)
        loc_system_.row_dofs[loc_ele_dof] = loc_system_.col_dofs[loc_ele_dof] = global_dofs_[loc_ele_dof];
        
        //shortcuts
        const unsigned int nsides = ele->n_sides();
        LinSys *ls = ad_->lin_sys;
        
        unsigned int side_row, edge_row;
        
        dirichlet_edge.resize(nsides);
        for (unsigned int i = 0; i < nsides; i++) {

            side_row = loc_side_dofs[i];    //local
            edge_row = loc_edge_dofs[i];    //local
            loc_system_.row_dofs[side_row] = loc_system_.col_dofs[side_row] = global_dofs_[side_row];    //global
            loc_system_.row_dofs[edge_row] = loc_system_.col_dofs[edge_row] = global_dofs_[edge_row];    //global
            
            dirichlet_edge[i] = 0;
            Side side = *dh_cell.elm().side(i);
            if (side.is_boundary()) {
                Boundary bcd = side.cond();
                ElementAccessor<3> b_ele = bcd.element_accessor();
                DarcyMH::EqData::BC_Type type = (DarcyMH::EqData::BC_Type)ad_->bc_type.value(b_ele.centre(), b_ele);

                double cross_section = ad_->cross_section.value(ele.centre(), ele);

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

                    unsigned int loc_edge_idx = bcd.bc_ele_idx();
                    char & switch_dirichlet = ad_->bc_switch_dirichlet[loc_edge_idx];
                    double bc_pressure = ad_->bc_switch_pressure.value(b_ele.centre(), b_ele);
                    double bc_flux = -ad_->bc_flux.value(b_ele.centre(), b_ele);
                    double side_flux = bc_flux * b_ele.measure() * cross_section;

                    // ** Update BC type. **
                    if (switch_dirichlet) {
                        // check and possibly switch to flux BC
                        // The switch raise error on the corresponding edge row.
                        // Magnitude of the error is abs(solution_flux - side_flux).
                        ASSERT_DBG(ad_->dh_->distr()->is_local(global_dofs_[side_row]))(global_dofs_[side_row]);
                        unsigned int loc_side_row = local_dofs_[side_row];
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
                        ASSERT_DBG(ad_->dh_->distr()->is_local(global_dofs_[edge_row]))(global_dofs_[edge_row]);
                        unsigned int loc_edge_row = local_dofs_[edge_row];
                        double & solution_head = ls->get_solution_array()[loc_edge_row];

                        if ( solution_head > bc_pressure) {
                            //DebugOut().fmt("x: {}, to dirich, p: {} -> p: {} f: {}\n",b_ele.centre()[0], solution_head, bc_pressure, bc_flux);
                            solution_head = bc_pressure;
                            switch_dirichlet=1;
                        }
                    }
                    
                        // ** Apply BCUpdate BC type. **
                        // Force Dirichlet type during the first iteration of the unsteady case.
                        if (switch_dirichlet || ad_->force_no_neumann_bc ) {
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
                    ASSERT_DBG(ad_->dh_->distr()->is_local(global_dofs_[edge_row]))(global_dofs_[edge_row]);
                    unsigned int loc_edge_row = local_dofs_[edge_row];
                    double & solution_head = ls->get_solution_array()[loc_edge_row];

                    // Force Robin type during the first iteration of the unsteady case.
                    if (solution_head > bc_switch_pressure  || ad_->force_no_neumann_bc) {
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
    }
        
     void assemble_sides(const DHCellAccessor& dh_cell)
    {
        const ElementAccessor<3> ele = dh_cell.elm();
        double cs = ad_->cross_section.value(ele.centre(), ele);
        double conduct =  ad_->conductivity.value(ele.centre(), ele);
        double scale = 1 / cs /conduct;
        
        assemble_sides_scale(dh_cell, scale);
    }
    
    void assemble_sides_scale(const DHCellAccessor& dh_cell, double scale)
    {
        arma::vec3 &gravity_vec = ad_->gravity_vec_;
        const ElementAccessor<3> ele = dh_cell.elm();
        
        fe_values_.reinit(ele);
        unsigned int ndofs = fe_values_.n_dofs();
        unsigned int qsize = fe_values_.n_points();
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
                                    (ad_->anisotropy.value(ele.centre(), ele )).i()
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
    
    
    void assemble_element(const DHCellAccessor&){
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
    

    void assembly_dim_connections(const DHCellAccessor& dh_cell){
        //D, E',E block: compatible connections: element-edge
        auto ele = dh_cell.elm(); //ElementAccessor<3>
        
        // no Neighbours => nothing to asssemble here
        if(dh_cell.elm()->n_neighs_vb() == 0) return;

        std::vector<LongIdx> higher_dim_dofs;
        //DebugOut() << "adc " << print_var(this) << print_var(side_quad_.size());
        unsigned int i = 0;
        for ( DHCellSide neighb_side : dh_cell.neighb_sides() ) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure

            loc_system_vb_.reset();
            loc_system_vb_.row_dofs[0] = loc_system_vb_.col_dofs[0] = global_dofs_[loc_ele_dof];
            
            // // read neighbor dofs (dh dofhandler)
            // // neighbor cell owning neighb_side
            DHCellAccessor dh_neighb_cell = neighb_side.cell();

            higher_dim_dofs.resize(dh_neighb_cell.n_dofs());
            dh_neighb_cell.get_dof_indices(higher_dim_dofs);

            // local index of pedge dof on neighboring cell
            // (dim+2) is number of edges of higher dim element
            // TODO: replace with DHCell getter when available for FESystem component
            const unsigned int t = dh_neighb_cell.n_dofs() - (dim+2) + neighb_side.side().side_idx();
            loc_system_vb_.row_dofs[1] = loc_system_vb_.col_dofs[1] = higher_dim_dofs[t];

            assembly_local_vb(ele, neighb_side);

            loc_system_vb_.eliminate_solution();
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

    void add_fluxes_in_balance_matrix(const DHCellAccessor& dh_cell){
        
        for(DHCellSide dh_side : dh_cell.side_range()){
            unsigned int sidx = dh_side.side_idx();
            if (dh_side.side().is_boundary()) {
                ad_->balance->add_flux_values(ad_->water_balance_idx, dh_side,
                                              {local_dofs_[loc_side_dofs[sidx]]},
                                              {1}, 0);
            }
        }
    }



    // assembly volume integrals
    FE_RT0<dim> fe_rt_;
    QGauss quad_;
    FEValues<3> fe_values_;

    NeighSideValues< (dim<3) ? dim : 2> ngh_values_;

    // Interpolation of velocity into barycenters
    QGauss velocity_interpolation_quad_;
    FEValues<3> velocity_interpolation_fv_;

    // data shared by assemblers of different dimension
    AssemblyDataPtrMH ad_;
    std::vector<unsigned int> dirichlet_edge;

    LocalSystem loc_system_;
    LocalSystem loc_system_vb_;
    std::vector<unsigned int> loc_side_dofs;
    std::vector<unsigned int> loc_edge_dofs;
    unsigned int loc_ele_dof;

    std::shared_ptr<MortarAssemblyBase> mortar_assembly;

    std::vector<LongIdx> global_dofs_;
    LocDofVec local_dofs_;
};


#endif /* SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_ */
