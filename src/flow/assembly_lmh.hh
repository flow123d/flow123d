/*
 * assembly_lmh.hh
 *
 *  Created on: Apr 21, 2016
 *      Author: jb
 */

#ifndef SRC_ASSEMBLY_LMH_HH_
#define SRC_ASSEMBLY_LMH_HH_

#include "system/index_types.hh"
#include "system/fmt/posix.h"
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

#include "la/linsys_PETSC.hh"
// #include "la/linsys_BDDC.hh"
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
template <int dim>
class AssemblyLMH : public AssemblyBase
{
public:
    typedef std::shared_ptr<DarcyLMH::EqData> AssemblyDataPtrLMH;
    
    AssemblyLMH<dim>(AssemblyDataPtrLMH data)
    : quad_(dim, 3),
      velocity_interpolation_quad_(dim, 0), // veloctiy values in barycenter
      ad_(data)
    {
        fe_values_.initialize(quad_, fe_rt_, update_values | update_gradients | update_JxW_values | update_quadrature_points);
        velocity_interpolation_fv_.initialize(velocity_interpolation_quad_, fe_rt_, update_values | update_quadrature_points);
        // local numbering of dofs for MH system
        // note: this shortcut supposes that the fe_system is the same on all elements
        // the function DiscreteSpace.fe(ElementAccessor) does not in fact depend on the element accessor
        auto fe = ad_->dh_->ds()->fe(ad_->dh_->own_range().begin()->elm())[Dim<dim>{}];
        FESystem<dim>* fe_system = dynamic_cast<FESystem<dim>*>(fe.get());
        loc_side_dofs = fe_system->fe_dofs(0);
        loc_ele_dof = fe_system->fe_dofs(1)[0];
        loc_edge_dofs = fe_system->fe_dofs(2);

        FEAL_ASSERT(ad_->mortar_method_ == DarcyFlowInterface::NoMortar)
            .error("Mortar methods are not supported in Lumped Mixed-Hybrid Method.");
        // if (ad_->mortar_method_ == DarcyFlowInterface::MortarP0) {
        //     mortar_assembly = std::make_shared<P0_CouplingAssembler>(ad_);
        // } else if (ad_->mortar_method_ == DarcyFlowInterface::MortarP1) {
        //     mortar_assembly = std::make_shared<P1_CouplingAssembler>(ad_);
        // }

        // reconstructed vector for the velocity and pressure
        // length = local Schur complement offset in LocalSystem
        schur_offset_ = loc_edge_dofs[0];
        reconstructed_solution_.zeros(schur_offset_);
    }


    ~AssemblyLMH<dim>() override
    {}
    
    void fix_velocity(const DHCellAccessor&) override
    {
        // if (mortar_assembly)
        //     mortar_assembly->fix_velocity(ele_ac);
    }

    void assemble_reconstruct(const DHCellAccessor& dh_cell) override
    {
        ASSERT_EQ_DBG(dh_cell.dim(), dim);
    
        assemble_local_system(dh_cell, false);   //do not switch dirichlet in seepage when reconstructing
        
        // TODO:
        // if (mortar_assembly)
        //     mortar_assembly->assembly(ele_ac);
        // if (mortar_assembly)
        //     mortar_assembly->fix_velocity(ele_ac);

        arma::vec schur_solution = ad_->p_edge_solution.get_subvec(loc_schur_.row_dofs);
        // reconstruct the velocity and pressure
        loc_system_.reconstruct_solution_schur(schur_offset_, schur_solution, reconstructed_solution_);

        // postprocess the velocity
        postprocess_velocity(dh_cell, reconstructed_solution_);
        
        ad_->full_solution.set_subvec(loc_system_.row_dofs.head(schur_offset_), reconstructed_solution_);
        ad_->full_solution.set_subvec(loc_system_.row_dofs.tail(loc_schur_.row_dofs.n_elem), schur_solution);
    }
    
    void assemble(const DHCellAccessor& dh_cell) override
    {
        ASSERT_EQ_DBG(dh_cell.dim(), dim);
        save_local_system_ = false;
        bc_fluxes_reconstruted = false;

        assemble_local_system(dh_cell, true);   //do use_dirichlet_switch
        
        loc_system_.compute_schur_complement(schur_offset_, loc_schur_, true);

        save_local_system(dh_cell);
        
        loc_schur_.eliminate_solution();
        ad_->lin_sys_schur->set_local_system(loc_schur_, ad_->dh_cr_->get_local_to_global_map());

        // TODO:
        // if (mortar_assembly)
        //     mortar_assembly->assembly(dh_cell);
    }

    void update_water_content(const DHCellAccessor&) override
    {};

protected:
    static unsigned int size()
    {
        // dofs: velocity, pressure, edge pressure
        return RefElement<dim>::n_sides + 1 + RefElement<dim>::n_sides;
    }
    
    void assemble_local_system(const DHCellAccessor& dh_cell, bool use_dirichlet_switch)
    {
        set_dofs(dh_cell);

        assemble_bc(dh_cell, use_dirichlet_switch);
        assemble_sides(dh_cell);
        assemble_element(dh_cell);
        assemble_source_term(dh_cell);
        assembly_dim_connections(dh_cell);
    }

    /** Loads the local system from a map: element index -> LocalSystem,
     * if it exits, or if the full solution is not yet reconstructed,
     * and reconstructs the full solution on the element.
     * Currently used only for seepage BC.
     */
    void load_local_system(const DHCellAccessor& dh_cell)
    {
        // do this only once per element
        if(bc_fluxes_reconstruted)
            return;

        // possibly find the corresponding local system
        auto ls = ad_->seepage_bc_systems.find(dh_cell.elm_idx());
        if (ls != ad_->seepage_bc_systems.end())
        {
            arma::vec schur_solution = ad_->p_edge_solution.get_subvec(loc_schur_.row_dofs);            
            // reconstruct the velocity and pressure
            ls->second.reconstruct_solution_schur(schur_offset_, schur_solution, reconstructed_solution_);

            postprocess_velocity(dh_cell, reconstructed_solution_);

            bc_fluxes_reconstruted = true;
        }
    }


    /// Saves the local system to a map: element index -> LocalSystem.
    /// Currently used only for seepage BC.
    void save_local_system(const DHCellAccessor& dh_cell)
    {
        // for seepage BC, save local system
        if(save_local_system_)
            ad_->seepage_bc_systems[dh_cell.elm_idx()] = loc_system_;
    }


    void set_dofs(const DHCellAccessor& dh_cell){
        const ElementAccessor<3> ele = dh_cell.elm();
        const DHCellAccessor dh_cr_cell = dh_cell.cell_with_other_dh(ad_->dh_cr_.get());

        unsigned int loc_size = size() + ele->n_neighs_vb();
        unsigned int loc_size_schur = ele->n_sides() + ele->n_neighs_vb();
        LocDofVec dofs(loc_size);
        LocDofVec dofs_schur(loc_size_schur);
        
        // add full vec indices
        dofs.head(dh_cell.n_dofs()) = dh_cell.get_loc_dof_indices();
        // add schur vec indices
        dofs_schur.head(dh_cr_cell.n_dofs()) = dh_cr_cell.get_loc_dof_indices();
        
        if(ele->n_neighs_vb() != 0)
        {
            //D, E',E block: compatible connections: element-edge
            unsigned int i = 0;
            
            for ( DHCellSide neighb_side : dh_cell.neighb_sides() ) {

                // read neighbor dofs (dh dofhandler)
                // neighbor cell owning neighb_side
                DHCellAccessor dh_neighb_cell = neighb_side.cell();
                
                // read neighbor dofs (dh_cr dofhandler)
                // neighbor cell owning neighb_side
                DHCellAccessor dh_neighb_cr_cell = dh_neighb_cell.cell_with_other_dh(ad_->dh_cr_.get());
                
                // local index of pedge dof in local system
                const unsigned int p = size()+i;
                // local index of pedge dof on neighboring cell
                const unsigned int t = dh_neighb_cell.n_dofs() - dh_neighb_cr_cell.n_dofs() + neighb_side.side().side_idx();
                dofs[p] = dh_neighb_cell.get_loc_dof_indices()[t];

                // local index of pedge dof in local schur system
                const unsigned int tt = dh_cr_cell.n_dofs()+i;
                dofs_schur[tt] = dh_neighb_cr_cell.get_loc_dof_indices()
                            [neighb_side.side().side_idx()];
                i++;
            }
        }
        loc_system_.reset(dofs, dofs);
        loc_schur_.reset(dofs_schur, dofs_schur);
    }

    
    void assemble_bc(const DHCellAccessor& dh_cell, bool use_dirichlet_switch){
        const ElementAccessor<3> ele = dh_cell.elm();
        
        dirichlet_edge.resize(ele->n_sides());
        for(DHCellSide dh_side : dh_cell.side_range()){
            unsigned int sidx = dh_side.side_idx();
            dirichlet_edge[sidx] = 0;

            // assemble BC
            if (dh_side.side().is_boundary()) {
                double cross_section = ad_->cross_section.value(ele.centre(), ele);
                assemble_side_bc(dh_side, cross_section, use_dirichlet_switch);

                ad_->balance->add_flux_values(ad_->water_balance_idx, dh_side,
                                              {loc_system_.row_dofs[loc_side_dofs[sidx]]},
                                              {1}, 0);
            }

            // side-edge (flux-lambda) terms
            loc_system_.add_value(loc_side_dofs[sidx], loc_edge_dofs[sidx], 1.0);
            loc_system_.add_value(loc_edge_dofs[sidx], loc_side_dofs[sidx], 1.0);
        }
    }


    void assemble_side_bc(const DHCellSide& side, double cross_section, bool use_dirichlet_switch)
    {
        const unsigned int sidx = side.side_idx();
        const unsigned int side_row = loc_side_dofs[sidx];    //local
        const unsigned int edge_row = loc_edge_dofs[sidx];    //local

        ElementAccessor<3> b_ele = side.cond().element_accessor();
        DarcyMH::EqData::BC_Type type = (DarcyMH::EqData::BC_Type)ad_->bc_type.value(b_ele.centre(), b_ele);

        if ( type == DarcyMH::EqData::none) {
            // homogeneous neumann
        } else if ( type == DarcyMH::EqData::dirichlet ) {
            double bc_pressure = ad_->bc_pressure.value(b_ele.centre(), b_ele);
            loc_schur_.set_solution(sidx, bc_pressure);
            dirichlet_edge[sidx] = 1;
            
        } else if ( type == DarcyMH::EqData::total_flux) {
            // internally we work with outward flux
            double bc_flux = -ad_->bc_flux.value(b_ele.centre(), b_ele);
            double bc_pressure = ad_->bc_pressure.value(b_ele.centre(), b_ele);
            double bc_sigma = ad_->bc_robin_sigma.value(b_ele.centre(), b_ele);
            
            dirichlet_edge[sidx] = 2;  // to be skipped in LMH source assembly
            loc_system_.add_value(edge_row, edge_row,
                                    -b_ele.measure() * bc_sigma * cross_section,
                                    (bc_flux - bc_sigma * bc_pressure) * b_ele.measure() * cross_section);
        }
        else if (type==DarcyMH::EqData::seepage) {
            ad_->is_linear=false;

            char & switch_dirichlet = ad_->bc_switch_dirichlet[b_ele.idx()];
            double bc_pressure = ad_->bc_switch_pressure.value(b_ele.centre(), b_ele);
            double bc_flux = -ad_->bc_flux.value(b_ele.centre(), b_ele);
            double side_flux = bc_flux * b_ele.measure() * cross_section;

            // ** Update BC type. **
            if(use_dirichlet_switch){  // skip BC change if only reconstructing the solution
            if (switch_dirichlet) {
                // check and possibly switch to flux BC
                // The switch raise error on the corresponding edge row.
                // Magnitude of the error is abs(solution_flux - side_flux).
                
                // try reconstructing local system for seepage BC
                load_local_system(side.cell());
                double solution_flux = reconstructed_solution_[side_row];

                if ( solution_flux < side_flux) {
                    //DebugOut().fmt("x: {}, to neum, p: {} f: {} -> f: {}\n", b_ele.centre()[0], bc_pressure, solution_flux, side_flux);
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

                double solution_head = ad_->p_edge_solution[loc_schur_.row_dofs[sidx]];

                if ( solution_head > bc_pressure) {
                    //DebugOut().fmt("x: {}, to dirich, p: {} -> p: {} f: {}\n",b_ele.centre()[0], solution_head, bc_pressure, bc_flux);
                    switch_dirichlet=1;
                }
            }
            }

            save_local_system_ = (bool) switch_dirichlet;
            
            // ** Apply BCUpdate BC type. **
            // Force Dirichlet type during the first iteration of the unsteady case.
            if (switch_dirichlet || ad_->force_no_neumann_bc ) {
                //DebugOut().fmt("x: {}, dirich, bcp: {}\n", b_ele.centre()[0], bc_pressure);
                loc_schur_.set_solution(sidx, bc_pressure);
                dirichlet_edge[sidx] = 1;
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

            double solution_head = ad_->p_edge_solution[loc_schur_.row_dofs[sidx]];

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


    virtual void assemble_sides(const DHCellAccessor& dh_cell)
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
        auto ele = dh_cell.elm();
        
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
                                    (ad_->anisotropy.value(ele.centre(), ele)).i()
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
//         if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
//             const arma::mat& local_matrix = loc_system_.get_matrix();
//             for(unsigned int i=0; i < ndofs; i++) {
//                 double val_side = local_matrix(i,i);
//                 double val_edge = -1./local_matrix(i,i);
// 
//                 unsigned int side_row = loc_system_.row_dofs[loc_side_dofs[i]];
//                 unsigned int edge_row = loc_system_.row_dofs[loc_edge_dofs[i]];
//                 static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( side_row, val_side );
//                 static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( edge_row, val_edge );
//             }
//         }
    }
    
    
    void assemble_element(FMT_UNUSED const DHCellAccessor& dh_cell){
        // set block B, B': element-side, side-element
        
        for(unsigned int side = 0; side < loc_side_dofs.size(); side++){
            loc_system_.add_value(loc_ele_dof, loc_side_dofs[side], -1.0);
            loc_system_.add_value(loc_side_dofs[side], loc_ele_dof, -1.0);
        }
        
//         if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
//             double val_ele =  1.;
//             static_cast<LinSys_BDDC*>(ad_->lin_sys)->
//                             diagonal_weights_set_value( loc_system_.row_dofs[loc_ele_dof], val_ele );
//         }
    }
    
    virtual void assemble_source_term(const DHCellAccessor& dh_cell)
    {
        const ElementAccessor<3> ele = dh_cell.elm();
        
        // compute lumped source
        double alpha = 1.0 / ele->n_sides();
        double cross_section = ad_->cross_section.value(ele.centre(), ele);
        double coef = alpha * ele.measure() * cross_section;
        
        double source = ad_->water_source_density.value(ele.centre(), ele)
                        + ad_->extra_source.value(ele.centre(), ele);
        double source_term = coef * source;
        
        // in unsteady, compute time term
        double storativity = 0.0;
        double time_term_diag = 0.0, time_term = 0.0, time_term_rhs = 0.0;
        
        if(! ad_->use_steady_assembly_)
        {
            storativity = ad_->storativity.value(ele.centre(), ele)
                          + ad_->extra_storativity.value(ele.centre(), ele);
            time_term = coef * storativity;
        }
        
        for (unsigned int i=0; i<ele->n_sides(); i++)
        {
            if(! ad_->use_steady_assembly_)
            {
                time_term_diag = time_term / ad_->time_step_;
                time_term_rhs = time_term_diag * ad_->p_edge_solution_previous_time[loc_schur_.row_dofs[i]];

                ad_->balance->add_mass_values(ad_->water_balance_idx, dh_cell,
                                              {loc_system_.row_dofs[loc_edge_dofs[i]]}, {time_term}, 0);
            }

            this->loc_system_.add_value(loc_edge_dofs[i], loc_edge_dofs[i],
                                        -time_term_diag,
                                        -source_term - time_term_rhs);

            ad_->balance->add_source_values(ad_->water_balance_idx, ele.region().bulk_idx(),
                                                {loc_system_.row_dofs[loc_edge_dofs[i]]}, {0},{source_term});
        }
    }

    void assembly_dim_connections(const DHCellAccessor& dh_cell){
        //D, E',E block: compatible connections: element-edge
        const ElementAccessor<3> ele = dh_cell.elm();
        
        // no Neighbours => nothing to asssemble here
        if(ele->n_neighs_vb() == 0) return;
        
        ASSERT_LT_DBG(ele->dim(), 3);
        arma::vec3 nv;

        unsigned int i = 0;
        for ( DHCellSide neighb_side : dh_cell.neighb_sides() ) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            unsigned int p = size()+i; // loc dof of higher ele edge
            
            ElementAccessor<3> ele_higher = neighb_side.cell().elm();
            ngh_values_.fe_side_values_.reinit(neighb_side.side());
            nv = ngh_values_.fe_side_values_.normal_vector(0);

            double value = ad_->sigma.value( ele.centre(), ele) *
                            2*ad_->conductivity.value( ele.centre(), ele) *
                            arma::dot(ad_->anisotropy.value( ele.centre(), ele)*nv, nv) *
                            ad_->cross_section.value( neighb_side.centre(), ele_higher ) * // cross-section of higher dim. (2d)
                            ad_->cross_section.value( neighb_side.centre(), ele_higher ) /
                            ad_->cross_section.value( ele.centre(), ele ) *      // crossection of lower dim.
                            neighb_side.measure();

            loc_system_.add_value(loc_ele_dof, loc_ele_dof, -value);
            loc_system_.add_value(loc_ele_dof, p,            value);
            loc_system_.add_value(p,loc_ele_dof,             value);
            loc_system_.add_value(p,p,                      -value);

//             // update matrix for weights in BDDCML
//             if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
//                 int ind = loc_system_.row_dofs[p];
//                // there is -value on diagonal in block C!
//                static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( ind, -value );
//             }
            ++i;
        }
    }


    void postprocess_velocity(const DHCellAccessor& dh_cell, arma::vec& solution)
    {
        const ElementAccessor<3> ele = dh_cell.elm();
        
        double edge_scale = ele.measure()
                              * ad_->cross_section.value(ele.centre(), ele)
                              / ele->n_sides();
        
        double edge_source_term = edge_scale * 
                ( ad_->water_source_density.value(ele.centre(), ele)
                + ad_->extra_source.value(ele.centre(), ele));
      
        postprocess_velocity_specific(dh_cell, solution, edge_scale, edge_source_term);
    }

    virtual void postprocess_velocity_specific(const DHCellAccessor& dh_cell, arma::vec& solution,
                                               double edge_scale, double edge_source_term)// override
    {
        const ElementAccessor<3> ele = dh_cell.elm();
        
        double storativity = ad_->storativity.value(ele.centre(), ele)
                             + ad_->extra_storativity.value(ele.centre(), ele);
        double new_pressure, old_pressure, time_term = 0.0;
        
        for (unsigned int i=0; i<ele->n_sides(); i++) {
            
            if( ! ad_->use_steady_assembly_)
            {
                new_pressure = ad_->p_edge_solution[loc_schur_.row_dofs[i]];
                old_pressure = ad_->p_edge_solution_previous_time[loc_schur_.row_dofs[i]];
                time_term = edge_scale * storativity / ad_->time_step_ * (new_pressure - old_pressure);
            }
            solution[loc_side_dofs[i]] += edge_source_term - time_term;
        }
    }

    // assembly volume integrals
    FE_RT0<dim> fe_rt_;
    QGauss quad_;
    FEValues<3> fe_values_;

    NeighSideValues<dim<3?dim:2> ngh_values_;

    // Interpolation of velocity into barycenters
    QGauss velocity_interpolation_quad_;
    FEValues<3> velocity_interpolation_fv_;

    // data shared by assemblers of different dimension
    AssemblyDataPtrLMH ad_;
    
    /** TODO: Investigate why the hell do we need this flag.
    *  If removed, it does not break any of the integration tests,
    * however it must influence the Dirichlet rows in matrix.
    */
    std::vector<unsigned int> dirichlet_edge;

    LocalSystem loc_system_;
    LocalSystem loc_schur_;
    std::vector<unsigned int> loc_side_dofs;
    std::vector<unsigned int> loc_edge_dofs;
    unsigned int loc_ele_dof;

    // std::shared_ptr<MortarAssemblyBase> mortar_assembly;

    /// Index offset in the local system for the Schur complement.
    unsigned int schur_offset_;

    /// Vector for reconstruted solution (velocity and pressure on element) from Schur complement.
    arma::vec reconstructed_solution_;

    /// Flag for saving the local system.
    /// Currently used only in case of seepage BC.
    bool save_local_system_;

    /// Flag indicating whether the fluxes for seepage BC has been reconstructed already.
    bool bc_fluxes_reconstruted;
};


#endif /* SRC_ASSEMBLY_LMH_HH_ */
