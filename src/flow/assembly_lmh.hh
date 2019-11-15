/*
 * assembly_lmh.hh
 *
 *  Created on: Apr 21, 2016
 *      Author: jb
 */

#ifndef SRC_ASSEMBLY_LMH_HH_
#define SRC_ASSEMBLY_LMH_HH_

#include "system/index_types.hh"
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
        fe_values_(map_, quad_, fe_rt_,
                update_values | update_gradients | update_JxW_values | update_quadrature_points),

        velocity_interpolation_quad_(dim, 0), // veloctiy values in barycenter
        velocity_interpolation_fv_(map_,velocity_interpolation_quad_, fe_rt_, update_values | update_quadrature_points),

        ad_(data),
        loc_system_(size(), size())
    {
        // local numbering of dofs for MH system
        // note: this shortcut supposes that the fe_system is the same on all elements
        // the function DiscreteSpace.fe(ElementAccessor) does not in fact depend on the element accessor
        auto fe = ad_->dh_->ds()->fe(ad_->dh_->own_range().begin()->elm()).get<dim>();
        FESystem<dim>* fe_system = dynamic_cast<FESystem<dim>*>(fe.get());
        loc_side_dofs = fe_system->fe_dofs(0);
        loc_ele_dof = fe_system->fe_dofs(1)[0];
        loc_edge_dofs = fe_system->fe_dofs(2);
        
        unsigned int nsides = dim+1;
        
        // create local sparsity pattern
        base_local_sp_.set_size(size(),size());
        base_local_sp_.zeros();
        base_local_sp_.submat(0, 0, nsides, nsides).ones();
        base_local_sp_.diag().ones();
        // armadillo 8.4.3 bug with negative sub diagonal index
        // sp.diag(nsides+1).ones();
        // sp.diag(-nsides-1).ones();
        // sp.print();
        
        base_local_sp_.submat(0, nsides+1, nsides-1, size()-1).diag().ones();
        base_local_sp_.submat(nsides+1, 0, size()-1, nsides-1).diag().ones();
        
        loc_system_.set_sparsity(base_local_sp_);
        
        loc_schur_.reset(nsides, nsides);
        schur_sp_.ones(nsides, nsides);
        loc_schur_.set_sparsity(schur_sp_);

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

//     LocalSystem& get_local_system() override
//         { return loc_system_;}
    
    void fix_velocity(LocalElementAccessorBase<3> ele_ac) override
    {
        // if (mortar_assembly)
        //     mortar_assembly->fix_velocity(ele_ac);
    }

    void assemble_reconstruct(LocalElementAccessorBase<3> ele_ac) override
    {
        ASSERT_EQ_DBG(ele_ac.dim(), dim);
        loc_system_.reset();
        reconstruct = true;
    
        LocDofVec dofs, dofs_schur;
        set_loc_dofs_vec(ele_ac, dofs, dofs_schur);
        loc_system_.reset(dofs.n_elem, dofs.n_elem);
        // set dofs for local schur, since these are used in assembly functions
        loc_schur_.reset(dofs_schur,dofs_schur);
        
        assemble_bc(ele_ac);
        
        assemble_sides(ele_ac);
        assemble_element(ele_ac);
        assemble_source_term(ele_ac);
        
        assembly_dim_connections(ele_ac);
        
        // TODO:
        // if (mortar_assembly)
        //     mortar_assembly->assembly(ele_ac);
        // if (mortar_assembly)
        //     mortar_assembly->fix_velocity(ele_ac);

        arma::vec schur_solution = ad_->p_edge_solution.get_subvec(dofs_schur);
        // reconstruct the velocity and pressure
        loc_system_.reconstruct_solution_schur(schur_offset_, schur_solution, reconstructed_solution_);

        // postprocess the velocity
        postprocess_velocity(ele_ac.dh_cell(), reconstructed_solution_);

        ad_->full_solution.set_subvec(dofs.head(schur_offset_), reconstructed_solution_);
        ad_->full_solution.set_subvec(dofs.tail(dofs_schur.n_elem), schur_solution);
    }
    
    void assemble(LocalElementAccessorBase<3> ele_ac) override
    {
        ASSERT_EQ_DBG(ele_ac.dim(), dim);
        reconstruct = false;
        save_local_system_ = false;
        bc_fluxes_reconstruted = false;

        set_dofs(ele_ac);
        assemble_bc(ele_ac);
        
        assemble_sides(ele_ac);
        assemble_element(ele_ac);
        assemble_source_term(ele_ac);
        
        assembly_dim_connections(ele_ac);
        
        loc_system_.compute_schur_complement(schur_offset_, loc_schur_, true);

        save_local_system(ele_ac.dh_cell());
        
        loc_schur_.eliminate_solution();
        ad_->lin_sys_schur->set_local_system(loc_schur_, ad_->dh_cr_->get_local_to_global_map());

        if (ad_->balance != nullptr)
            add_fluxes_in_balance_matrix(ele_ac);

        // TODO:
        // if (mortar_assembly)
        //     mortar_assembly->assembly(ele_ac);
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

    arma::vec3 make_element_vector(LocalElementAccessorBase<3> ele_ac) override
    {
        //START_TIMER("Assembly<dim>::make_element_vector");
        arma::vec3 flux_in_center;
        flux_in_center.zeros();
        auto ele = ele_ac.element_accessor();

        velocity_interpolation_fv_.reinit(ele);
        for (unsigned int li = 0; li < ele->n_sides(); li++) {
            flux_in_center += ad_->full_solution[ ele_ac.side_local_row(li) ]
                        * velocity_interpolation_fv_.vector_view(0).value(li,0);
        }


        flux_in_center /= ad_->cross_section.value(ele.centre(), ele );
        return flux_in_center;
    }

    void update_water_content(const DHCellAccessor& dh_cell) override
    {};

protected:
    static const unsigned int size()
    {
        // dofs: velocity, pressure, edge pressure
        return RefElement<dim>::n_sides + 1 + RefElement<dim>::n_sides;
    }

    void set_loc_dofs_vec(LocalElementAccessorBase<3> ele_ac, LocDofVec& dofs, LocDofVec& dofs_schur){
        ElementAccessor<3> ele = ele_ac.element_accessor();
        DHCellAccessor dh_cell = ele_ac.dh_cell();
        DHCellAccessor dh_cr_cell = dh_cell.cell_with_other_dh(ad_->dh_cr_.get());
        
        unsigned int loc_size = size() + ele->n_neighs_vb();
        unsigned int loc_size_schur = ele->n_sides() + ele->n_neighs_vb();
        // vector of DoFs
        dofs.set_size(loc_size);
        dofs_schur.set_size(loc_size_schur);
        
        // add full vec indices
        dofs.head(dh_cell.n_dofs()) = dh_cell.get_loc_dof_indices();
        // add schur vec indices
        dofs_schur.head(dh_cr_cell.n_dofs()) = dh_cr_cell.get_loc_dof_indices();
        
        if(ele->n_neighs_vb() == 0)
            return;

        //D, E',E block: compatible connections: element-edge
        Neighbour *ngh;
        unsigned int i = 0;
        
        for ( DHCellSide neighb_side : dh_cell.neighb_sides() ) {
            
            ngh = ele->neigh_vb[i];
            
            // read neighbor dofs (dh dofhandler)
            // neighbor cell owning neighb_side
            DHCellAccessor dh_neighb_cell = neighb_side.cell();
            
            // read neighbor dofs (dh_cr dofhandler)
            // neighbor cell owning neighb_side
            DHCellAccessor dh_neighb_cr_cell = dh_neighb_cell.cell_with_other_dh(ad_->dh_cr_.get());
            
            // find edge dofs of neighbor corresponding to neighb_side
            for (unsigned int j = 0; j < neighb_side.element().dim()+1; j++){
            	if (neighb_side.element()->edge_idx(j) == ngh->edge_idx()) {
                    unsigned int p = size()+i;
                    dofs[p] = dh_neighb_cell.get_loc_dof_indices()
                                [dh_neighb_cell.n_dofs() - dh_neighb_cr_cell.n_dofs() + j];
                    
                    unsigned int t = dh_cr_cell.n_dofs()+i;
                    dofs_schur[t] = dh_neighb_cr_cell.get_loc_dof_indices()[j];
                    break;
                }
            }
            i++;
        }
    }
    
    void set_dofs(LocalElementAccessorBase<3> ele_ac){
        ElementAccessor<3> ele = ele_ac.element_accessor();
        
        unsigned int loc_size = size() + ele->n_neighs_vb();
        unsigned int loc_size_schur = ele->n_sides() + ele->n_neighs_vb();
        LocDofVec dofs_schur(loc_size_schur);
        
        DHCellAccessor dh_cr_cell = ele_ac.dh_cell().cell_with_other_dh(ad_->dh_cr_.get());
        // add schur vec indices
        dofs_schur.head(dh_cr_cell.n_dofs()) = dh_cr_cell.get_loc_dof_indices();
        
        schur_sp_.ones(loc_size_schur, loc_size_schur);

        // if neighbor communication, then resize the local system and set dofs and sp
        local_sp_.zeros(loc_size, loc_size);
        local_sp_( 0,0, arma::size(base_local_sp_)) = base_local_sp_;
        
        
        if(ele->n_neighs_vb() != 0)
        {
            //D, E',E block: compatible connections: element-edge
            Neighbour *ngh;

            for (unsigned int i = 0; i < ele->n_neighs_vb(); i++) {
                // every compatible connection adds a 2x2 matrix involving
                // current element pressure  and a connected edge pressure
                ngh = ele_ac.element_accessor()->neigh_vb[i];
                LocalElementAccessorBase<3> acc_higher_dim( ele_ac.dh_cell().dh()->cell_accessor_from_element(ngh->edge()->side(0)->element().idx()) );
                
                DHCellAccessor higher_dh_cr_cell = acc_higher_dim.dh_cell().cell_with_other_dh(ad_->dh_cr_.get());
                
                for (unsigned int j = 0; j < ngh->edge()->side(0)->element().dim()+1; j++)
                    if (ngh->edge()->side(0)->element()->edge_idx(j) == ngh->edge_idx()) {
                        unsigned int p = size()+i;
                        // dofs[p] = acc_higher_dim.edge_row(j);
                        local_sp_(loc_ele_dof, p) = 1;
                        local_sp_(p, loc_ele_dof) = 1;
                        local_sp_(p, p) = 1;
                        
                        unsigned int t = dh_cr_cell.n_dofs()+i;
                        dofs_schur[t] = higher_dh_cr_cell.get_loc_dof_indices()[j];
                        break;
                    }
            }
        }
        loc_system_.reset(loc_size, loc_size);
        loc_system_.set_sparsity(local_sp_);
        
        loc_schur_.reset(dofs_schur, dofs_schur);
        loc_schur_.set_sparsity(schur_sp_);
    }
    
    void assemble_bc(LocalElementAccessorBase<3> ele_ac){
        //shortcuts
        const unsigned int nsides = ele_ac.n_sides();
        
        Boundary *bcd;
        unsigned int side_row, edge_row;
        
        dirichlet_edge.resize(nsides);
        for (unsigned int i = 0; i < nsides; i++) {

            side_row = loc_side_dofs[i];    //local
            edge_row = loc_edge_dofs[i];    //local
            
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
                    loc_schur_.set_solution(i, bc_pressure);
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
                    if(! reconstruct){  // skip BC change if only reconstructing the solution
                    if (switch_dirichlet) {
                        // check and possibly switch to flux BC
                        // The switch raise error on the corresponding edge row.
                        // Magnitude of the error is abs(solution_flux - side_flux).
                        
                        // try reconstructing local system for seepage BC
                        load_local_system(ele_ac.dh_cell());
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

                        double solution_head = ad_->p_edge_solution[loc_schur_.row_dofs[i]];

                        if ( solution_head > bc_pressure) {
                            //DebugOut().fmt("x: {}, to dirich, p: {} -> p: {} f: {}\n",b_ele.centre()[0], solution_head, bc_pressure, bc_flux);
                            switch_dirichlet=1;
                        }
                    }
                    }

                    save_local_system_ = (bool) switch_dirichlet;
                    
                    // ** Apply BCUpdate BC type. **
                    // Force Dirichlet type during the first iteration of the unsteady case.
                    if (switch_dirichlet || ad_->force_bc_switch ) {
                        //DebugOut().fmt("x: {}, dirich, bcp: {}\n", b_ele.centre()[0], bc_pressure);
                        loc_schur_.set_solution(i, bc_pressure);
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

                    double solution_head = ad_->p_edge_solution[loc_schur_.row_dofs[i]];

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
    }
        
    void assemble_sides(LocalElementAccessorBase<3> ele_ac) override
    {
        ElementAccessor<3> ele = ele_ac.element_accessor();
        double cs = ad_->cross_section.value(ele.centre(), ele);
        double conduct =  ad_->conductivity.value(ele.centre(), ele);
        double scale = 1 / cs /conduct;
        
        assemble_sides_scale(ele_ac, scale);
    }
    
    void assemble_sides_scale(LocalElementAccessorBase<3> ele_ac, double scale)
    {
        arma::vec3 &gravity_vec = ad_->gravity_vec_;
        
        ElementAccessor<3> ele = ele_ac.element_accessor();
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
    
    
    void assemble_element(LocalElementAccessorBase<3> ele_ac){
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
    
    void assemble_source_term(LocalElementAccessorBase<3> ele_ac) override
    {
        ElementAccessor<3> ele = ele_ac.element_accessor();
        
        // compute lumped source
        double alpha = 1.0 / ele->n_sides();
        double cross_section = ad_->cross_section.value(ele.centre(), ele);
        double coef = alpha * ele.measure() * cross_section;
        
        double source = ad_->water_source_density.value(ele.centre(), ele);
        double source_term = coef * source;
        
        // in unsteady, compute time term
        double storativity = 0.0;
        double time_term_diag = 0.0, time_term = 0.0, time_term_rhs = 0.0;
        
        if(! ad_->use_steady_assembly_)
        {
            storativity = ad_->storativity.value(ele.centre(), ele);
            time_term = coef * storativity;
        }
        
        for (unsigned int i=0; i<ele->n_sides(); i++)
        {
            if(! ad_->use_steady_assembly_)
            {
                time_term_diag = time_term / ad_->time_step_;
                time_term_rhs = time_term_diag * ad_->p_edge_solution_previous_time[loc_schur_.row_dofs[i]];
            }

            this->loc_system_.add_value(loc_edge_dofs[i], loc_edge_dofs[i],
                                        -time_term_diag,
                                        -source_term - time_term_rhs);
            
            if( ! reconstruct)
            if (ad_->balance != nullptr)
            {
                ad_->balance->add_source_values(ad_->water_balance_idx, ele.region().bulk_idx(),
                                                {(LongIdx)ele_ac.edge_local_row(i)}, {0},{source_term});
                if( ! ad_->use_steady_assembly_)
                {
                    ad_->balance->add_mass_matrix_values(ad_->water_balance_idx, ele.region().bulk_idx(),
                                                {(LongIdx)ele_ac.edge_row(i)}, {time_term});
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
        
        ASSERT_LT_DBG(ele->dim(), 3);
        arma::vec3 nv;

        unsigned int i = 0;
        for ( DHCellSide neighb_side : dh_cell.neighb_sides() ) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            unsigned int p = size()+i; // loc dof of higher ele edge
            
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

    void add_fluxes_in_balance_matrix(LocalElementAccessorBase<3> ele_ac){

        ElementAccessor<3> ele = ele_ac.element_accessor();
        
        for (unsigned int i = 0; i < ele->n_sides(); i++) {
            Boundary* bcd = ele.side(i)->cond();

            if (bcd) {
                ad_->balance->add_flux_matrix_values(ad_->water_balance_idx, ele_ac.side(i),
                                                     {(LongIdx)(ele_ac.side_row(i))}, {1});
            }
        }
    }

    void postprocess_velocity(const DHCellAccessor& dh_cell, arma::vec& solution)
    {
        ElementAccessor<3> ele = dh_cell.elm();
        
        double edge_scale = ele.measure()
                              * ad_->cross_section.value(ele.centre(), ele)
                              / ele->n_sides();
        
        double edge_source_term = edge_scale * ad_->water_source_density.value(ele.centre(), ele);
      
        postprocess_velocity_specific(dh_cell, solution, edge_scale, edge_source_term);
    }

    virtual void postprocess_velocity_specific(const DHCellAccessor& dh_cell, arma::vec& solution,
                                               double edge_scale, double edge_source_term)// override
    {
        ElementAccessor<3> ele = dh_cell.elm();
        
        double storativity = ad_->storativity.value(ele.centre(), ele);
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
    
    
    // temporary fix in schur reconstruction
    bool reconstruct;

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
    arma::umat base_local_sp_;      ///< Sparsity pattern of the LocalSystem (without dim communication).
    arma::umat local_sp_;           ///< Whole sparsity pattern of the LocalSystem.
    arma::umat schur_sp_;           ///< Whole sparsity pattern of the Schur complement of LocalSystem.
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
