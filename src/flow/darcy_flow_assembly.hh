/*
 * darcy_flow_assembly.hh
 *
 *  Created on: Apr 21, 2016
 *      Author: jb
 */

#ifndef SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_
#define SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_

#include "mesh/long_idx.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/neighbours.h"
#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_values_views.hh"
#include "quadrature/quadrature_lib.hh"
#include "flow/mh_dofhandler.hh"

#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "la/linsys_BDDC.hh"
#include "la/schur.hh"

#include "la/local_system.hh"

#include "coupling/balance.hh"
#include "flow/mortar_assembly.hh"


class AssemblyBase
{
public:
    typedef std::shared_ptr<DarcyMH::EqData> AssemblyDataPtr;
    typedef std::vector<std::shared_ptr<AssemblyBase> > MultidimAssembly;

    virtual ~AssemblyBase() {}

    /**
        * Generic creator of multidimensional assembly, i.e. vector of
        * particular assembly objects.
        */
    template< template<int dim> class Impl >
    static MultidimAssembly create(typename Impl<1>::AssemblyDataPtr data) {
        return { std::make_shared<Impl<1> >(data),
            std::make_shared<Impl<2> >(data),
            std::make_shared<Impl<3> >(data) };

    }

//     virtual LocalSystem & get_local_system() = 0;
    virtual void fix_velocity(LocalElementAccessorBase<3> ele_ac) = 0;
    virtual void assemble(LocalElementAccessorBase<3> ele_ac) = 0;
    virtual void assemble_reconstruct(LocalElementAccessorBase<3> ele_ac) = 0;
        
    // assembly compatible neighbourings
//     virtual void assembly_local_vb(ElementAccessor<3> ele, Neighbour *ngh) = 0;

    // compute velocity value in the barycenter
    // TODO: implement and use general interpolations between discrete spaces
    virtual arma::vec3 make_element_vector(LocalElementAccessorBase<3> ele_ac) = 0;

    virtual void update_water_content(LocalElementAccessorBase<3> ele)
    {}

protected:

    virtual void assemble_sides(LocalElementAccessorBase<3> ele) =0;
    
    virtual void assemble_source_term(LocalElementAccessorBase<3> ele)
    {}
};



template <int dim>
class NeighSideValues {
private:
    // assembly face integrals (BC)
    MappingP1<dim+1,3> side_map_;
    QGauss<dim> side_quad_;
    FE_P_disc<dim+1> fe_p_disc_;
public:
    NeighSideValues<dim>()
    :  side_quad_(1),
       fe_p_disc_(0),
       fe_side_values_(side_map_, side_quad_, fe_p_disc_, update_normal_vectors)
    {}
    FESideValues<dim+1,3> fe_side_values_;

};



template<int dim>
class AssemblyMH : public AssemblyBase
{
public:
    AssemblyMH<dim>(AssemblyDataPtr data)
    : quad_(3),
        fe_values_(map_, quad_, fe_rt_,
                update_values | update_gradients | update_JxW_values | update_quadrature_points),

        velocity_interpolation_quad_(0), // veloctiy values in barycenter
        velocity_interpolation_fv_(map_,velocity_interpolation_quad_, fe_rt_, update_values | update_quadrature_points),

        ad_(data),
        loc_system_(size(), size())
    {
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

        if (ad_->mortar_method_ == DarcyMH::MortarP0) {
            mortar_assembly = std::make_shared<P0_CouplingAssembler>(ad_);
        } else if (ad_->mortar_method_ == DarcyMH::MortarP1) {
            mortar_assembly = std::make_shared<P1_CouplingAssembler>(ad_);
        }

    }


    ~AssemblyMH<dim>() override
    {}

//     LocalSystem& get_local_system() override
//         { return loc_system_;}
    
    void fix_velocity(LocalElementAccessorBase<3> ele_ac) override
    {
        if (mortar_assembly)
            mortar_assembly->fix_velocity(ele_ac);
    }

    void assemble_reconstruct(LocalElementAccessorBase<3> ele_ac) override
    {
        ASSERT_EQ_DBG(ele_ac.dim(), dim);
        loc_system_.reset();
        compute_balance = false;
    
        arma::Col<LongIdx> dofs, dofs_schur;
//         set_loc_dofs_vec(ele_ac, loc_system_.row_dofs, loc_schur_.row_dofs);
        set_loc_dofs_vec(ele_ac, dofs, dofs_schur);
        loc_system_.reset(dofs, dofs);
        loc_schur_.reset(dofs_schur,dofs_schur);
        
        assemble_bc(ele_ac);
        
        assemble_sides(ele_ac);
        assemble_element(ele_ac);
        assemble_source_term(ele_ac);
        
        assembly_dim_connections(ele_ac);
        
//         loc_system_.eliminate_solution();
        
        const unsigned int n = loc_schur_.row_dofs.n_rows;  // local Schur complement size
        const unsigned int off = loc_edge_dofs[0];          // local Schur complement offset in LocalSystem
        arma::vec schur_solution(n);
        for(unsigned int i=0; i<n; i++)
        {
            // read the computed edge pressures
            schur_solution(i) = (*ad_->schur_solution_)[dofs_schur[i]];
            // write the computed edge pressures
            unsigned int loc_row = dofs[off+i];
            (*ad_->full_solution_)[loc_row] = schur_solution(i);
        }
        
        // reconstruct the velocity and pressure
        arma::vec reconstructed_solution(off);
        loc_system_.reconstruct_solution_schur(off, schur_solution, reconstructed_solution);
        
        for(unsigned int i=0; i<off; i++)
        {
            // write the velocity and pressure
            unsigned int loc_row = dofs[i];
            (*ad_->full_solution_)[loc_row] += reconstructed_solution(i);
        }
        
//         if (mortar_assembly)
//             mortar_assembly->assembly(ele_ac);
    }
    
    void assemble(LocalElementAccessorBase<3> ele_ac) override
    {
        ASSERT_EQ_DBG(ele_ac.dim(), dim);
        compute_balance = true;
    
        set_dofs(ele_ac);
        assemble_bc(ele_ac);
        
        assemble_sides(ele_ac);
        assemble_element(ele_ac);
        assemble_source_term(ele_ac);
        
        assembly_dim_connections(ele_ac);
        
        ad_->lin_sys->set_local_system(loc_system_);
        
//         local.eliminate_solution();
        loc_system_.compute_schur_complement(loc_edge_dofs[0], loc_schur_, true);
        ad_->lin_sys_schur->set_local_system(loc_schur_);

        if (ad_->balance != nullptr)
            add_fluxes_in_balance_matrix(ele_ac);

        if (mortar_assembly)
            mortar_assembly->assembly(ele_ac);
    }

//     void assembly_local_vb(ElementAccessor<3> ele, Neighbour *ngh) override
//     {
//         ASSERT_LT_DBG(ele->dim(), 3);
//         //DebugOut() << "alv " << print_var(this);
//         //START_TIMER("Assembly<dim>::assembly_local_vb");
//         // compute normal vector to side
//         arma::vec3 nv;
//         ElementAccessor<3> ele_higher = ad_->mesh->element_accessor( ngh->side()->element().idx() );
//         ngh_values_.fe_side_values_.reinit(ele_higher, ngh->side()->side_idx());
//         nv = ngh_values_.fe_side_values_.normal_vector(0);
// 
//         double value = ad_->sigma.value( ele.centre(), ele) *
//                         2*ad_->conductivity.value( ele.centre(), ele) *
//                         arma::dot(ad_->anisotropy.value( ele.centre(), ele)*nv, nv) *
//                         ad_->cross_section.value( ngh->side()->centre(), ele_higher ) * // cross-section of higher dim. (2d)
//                         ad_->cross_section.value( ngh->side()->centre(), ele_higher ) /
//                         ad_->cross_section.value( ele.centre(), ele ) *      // crossection of lower dim.
//                         ngh->side()->measure();
// 
//         loc_system_vb_.add_value(0,0, -value);
//         loc_system_vb_.add_value(0,1,  value);
//         loc_system_vb_.add_value(1,0,  value);
//         loc_system_vb_.add_value(1,1, -value);
//     }


    arma::vec3 make_element_vector(LocalElementAccessorBase<3> ele_ac) override
    {
        //START_TIMER("Assembly<dim>::make_element_vector");
        arma::vec3 flux_in_center;
        flux_in_center.zeros();
        auto ele = ele_ac.element_accessor();

        velocity_interpolation_fv_.reinit(ele);
        for (unsigned int li = 0; li < ele->n_sides(); li++) {
            flux_in_center += ad_->data_vec_[ ele_ac.side_row(li) ] //ad_->mh_dh->side_flux( *(ele.side( li ) ) )
                        * velocity_interpolation_fv_.vector_view(0).value(li,0);
        }


        flux_in_center /= ad_->cross_section.value(ele.centre(), ele );
        return flux_in_center;
    }

protected:
    static const unsigned int size()
    {
        // dofs: velocity, pressure, edge pressure
        return RefElement<dim>::n_sides + 1 + RefElement<dim>::n_sides;
    }
    
    void set_loc_dofs_vec(LocalElementAccessorBase<3> ele_ac, arma::Col<LongIdx>& dofs, arma::Col<LongIdx>& dofs_schur){
        ElementAccessor<3> ele = ele_ac.element_accessor();
        
        unsigned int loc_size = size() + ele->n_neighs_vb();
        unsigned int loc_size_schur = ele_ac.n_sides() + ele->n_neighs_vb();
        // vector of DoFs
        dofs.set_size(loc_size);
        dofs_schur.set_size(loc_size_schur);
        
        std::vector<LongIdx> indices(ele_ac.dh_cell().n_dofs());
        ele_ac.dh_cell().get_loc_dof_indices(indices);
        
        DHCellAccessor dh_cr_cell = ele_ac.dh_cell().cell_with_other_dh(ad_->dh_cr_.get());
        std::vector<LongIdx> schur_indices(dh_cr_cell.n_dofs());
        dh_cr_cell.get_loc_dof_indices(schur_indices);
        
        if(ele->n_neighs_vb() == 0)
        {   
            dofs = indices;
            dofs_schur = schur_indices;
            return;
        }
        
        // add full vec indices
        arma::Col<LongIdx> d = indices;
        dofs.head(indices.size()) = d;
        
        // add schur vec indices
        arma::Col<LongIdx> dd = schur_indices;
        dofs_schur.head(schur_indices.size()) = dd;
        
        //D, E',E block: compatible connections: element-edge
        Neighbour *ngh;

        for (unsigned int i = 0; i < ele->n_neighs_vb(); i++) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            ngh = ele_ac.element_accessor()->neigh_vb[i];
            LocalElementAccessorBase<3> acc_higher_dim( ele_ac.dh_cell().dh()->cell_accessor_from_element(ngh->edge()->side(0)->element().idx()) );
            
            std::vector<LongIdx> higher_indices(acc_higher_dim.dh_cell().n_dofs());
            acc_higher_dim.dh_cell().get_loc_dof_indices(higher_indices);
            
            DHCellAccessor higher_dh_cr_cell = acc_higher_dim.dh_cell().cell_with_other_dh(ad_->dh_cr_.get());
            std::vector<LongIdx> cr_dofs(higher_dh_cr_cell.n_dofs());
            higher_dh_cr_cell.get_loc_dof_indices(cr_dofs);            
                    
            for (unsigned int j = 0; j < ngh->edge()->side(0)->element().dim()+1; j++)
            	if (ngh->edge()->side(0)->element()->edge_idx(j) == ngh->edge_idx()) {
                    unsigned int p = size()+i;
                    dofs[p] = higher_indices[higher_indices.size() - cr_dofs.size() + j];
                    
                    unsigned int t = schur_indices.size()+i;
                    dofs_schur[t] = cr_dofs[j];
            		break;
            	}
        }
    }
    
    void set_dofs(LocalElementAccessorBase<3> ele_ac){
        ElementAccessor<3> ele = ele_ac.element_accessor();
        
        unsigned int loc_size = size() + ele->n_neighs_vb();
        unsigned int loc_size_schur = ele_ac.n_sides() + ele->n_neighs_vb();
        // vector of DoFs
        arma::Col<LongIdx> dofs(loc_size);
        arma::Col<LongIdx> dofs_schur(loc_size_schur);
        
        DHCellAccessor dh_cr_cell = ele_ac.dh_cell().cell_with_other_dh(ad_->dh_cr_.get());
        std::vector<LongIdx> indices(dh_cr_cell.n_dofs());
        dh_cr_cell.get_dof_indices(indices);
        
        //set global dof for element (pressure)
        dofs[loc_ele_dof] = ele_ac.ele_row();
        
        unsigned int side_row, edge_row;
        for (unsigned int i = 0; i < ele_ac.n_sides(); i++) {

            side_row = loc_side_dofs[i];    //local
            edge_row = loc_edge_dofs[i];    //local
            
            dofs[side_row] = ele_ac.side_row(i);    //global
            dofs[edge_row] = ele_ac.edge_row(i);    //global
        }
        
        if(ele->n_neighs_vb() == 0)
        {
            loc_system_.reset(dofs, dofs);
            loc_system_.set_sparsity(base_local_sp_);
            
            dofs_schur = indices;
            loc_schur_.reset(dofs_schur,dofs_schur);
            schur_sp_.ones(loc_size_schur, loc_size_schur);
            loc_schur_.set_sparsity(schur_sp_);
            return;
        }
        
        // if neighbor communication, then resize the local system and set dofs and sp
        local_sp_.zeros(loc_size, loc_size);
        local_sp_( 0,0, arma::size(base_local_sp_)) = base_local_sp_;
        
        // TODO: sparse for commmunication blocks?
        schur_sp_.ones(loc_size_schur, loc_size_schur);
        
       
        arma::Col<LongIdx> d = indices;
        dofs_schur.head(indices.size()) = d;
        
        //D, E',E block: compatible connections: element-edge
        Neighbour *ngh;

        for (unsigned int i = 0; i < ele->n_neighs_vb(); i++) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            ngh = ele_ac.element_accessor()->neigh_vb[i];
            LocalElementAccessorBase<3> acc_higher_dim( ele_ac.dh_cell().dh()->cell_accessor_from_element(ngh->edge()->side(0)->element().idx()) );
            
            DHCellAccessor higher_dh_cr_cell = acc_higher_dim.dh_cell().cell_with_other_dh(ad_->dh_cr_.get());
            std::vector<int> cr_dofs(higher_dh_cr_cell.n_dofs());
            higher_dh_cr_cell.get_dof_indices(cr_dofs);
                    
            for (unsigned int j = 0; j < ngh->edge()->side(0)->element().dim()+1; j++)
            	if (ngh->edge()->side(0)->element()->edge_idx(j) == ngh->edge_idx()) {
                    unsigned int p = size()+i;
                    dofs[p] = acc_higher_dim.edge_row(j);
                    local_sp_(loc_ele_dof, p) = 1;
                    local_sp_(p, loc_ele_dof) = 1;
                    local_sp_(p, p) = 1;
                    
                    unsigned int t = indices.size()+i;
                    dofs_schur[t] = cr_dofs[j];
            		break;
            	}
        }
        
        loc_system_.reset(dofs, dofs);
        loc_system_.set_sparsity(local_sp_);
        
        loc_schur_.reset(dofs_schur, dofs_schur);
        loc_schur_.set_sparsity(schur_sp_);
    }
    
    void assemble_bc(LocalElementAccessorBase<3> ele_ac){
        //shortcuts
        const unsigned int nsides = ele_ac.n_sides();
        LinSys *ls = ad_->lin_sys;
        
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
    
    void assemble_source_term(LocalElementAccessorBase<3> ele_ac)
    {
        double cs = ad_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor());

        // set sources
        double source = ele_ac.measure() * cs *
                ad_->water_source_density.value(ele_ac.centre(), ele_ac.element_accessor());
        loc_system_.add_value(loc_ele_dof, -1.0 * source);

        if(compute_balance)
            ad_->balance->add_source_values(ad_->water_balance_idx, ele_ac.region().bulk_idx(), {(LongIdx) ele_ac.ele_local_idx()}, {0}, {source});
    }

    void assembly_dim_connections(LocalElementAccessorBase<3> ele_ac){
        //D, E',E block: compatible connections: element-edge
        ElementAccessor<3> ele = ele_ac.element_accessor();
        
        // no Neighbours => nothing to asssemble here
        if(ele->n_neighs_vb() == 0) return;
        
        ASSERT_LT_DBG(ele->dim(), 3);
        int ele_row = ele_ac.ele_row();
        Neighbour *ngh;
        arma::vec3 nv;

        //DebugOut() << "adc " << print_var(this) << print_var(side_quad_.size());
        for (unsigned int i = 0; i < ele->n_neighs_vb(); i++) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            unsigned int p = size()+i; // loc dof of higher ele edge
            ngh = ele_ac.element_accessor()->neigh_vb[i];
            
            //DebugOut() << "alv " << print_var(this);
            // compute normal vector to side
            ElementAccessor<3> ele_higher = ad_->mesh->element_accessor( ngh->side()->element().idx() );
            ngh_values_.fe_side_values_.reinit(ele_higher, ngh->side()->side_idx());
            nv = ngh_values_.fe_side_values_.normal_vector(0);

            double value = ad_->sigma.value( ele.centre(), ele) *
                            2*ad_->conductivity.value( ele.centre(), ele) *
                            arma::dot(ad_->anisotropy.value( ele.centre(), ele)*nv, nv) *
                            ad_->cross_section.value( ngh->side()->centre(), ele_higher ) * // cross-section of higher dim. (2d)
                            ad_->cross_section.value( ngh->side()->centre(), ele_higher ) /
                            ad_->cross_section.value( ele.centre(), ele ) *      // crossection of lower dim.
                            ngh->side()->measure();

            loc_system_.add_value(loc_ele_dof, loc_ele_dof, -value);
            loc_system_.add_value(loc_ele_dof, p,            value);
            loc_system_.add_value(p,loc_ele_dof,             value);
            loc_system_.add_value(p,p,                      -value);

            // update matrix for weights in BDDCML
            if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
                int ind = loc_system_.row_dofs[p];
               // there is -value on diagonal in block C!
               static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( ind, -value );
            }
        }
    }

    void add_fluxes_in_balance_matrix(LocalElementAccessorBase<3> ele_ac){

        for (unsigned int i = 0; i < ele_ac.n_sides(); i++) {
            Boundary* bcd = ele_ac.side(i)->cond();

            if (bcd) {
                ad_->balance->add_flux_matrix_values(ad_->water_balance_idx, ad_->local_boundary_index,
                                                     {(LongIdx)(ele_ac.side_row(i))}, {1});
                ++(ad_->local_boundary_index);
            }
        }
    }


    // temporary fix in schur reconstruction
    bool compute_balance;

    // assembly volume integrals
    FE_RT0<dim> fe_rt_;
    MappingP1<dim,3> map_;
    QGauss<dim> quad_;
    FEValues<dim,3> fe_values_;

    NeighSideValues<dim<3?dim:2> ngh_values_;

    // Interpolation of velocity into barycenters
    QGauss<dim> velocity_interpolation_quad_;
    FEValues<dim,3> velocity_interpolation_fv_;

    // data shared by assemblers of different dimension
    AssemblyDataPtr ad_;
    std::vector<unsigned int> dirichlet_edge;

    LocalSystem loc_system_;
    arma::umat base_local_sp_;      ///< Sparsity pattern of the LocalSystem (without dim communication).
    arma::umat local_sp_;           ///< Whole sparsity pattern of the LocalSystem.
    arma::umat schur_sp_;           ///< Whole sparsity pattern of the Schur complement of LocalSystem.
    LocalSystem loc_schur_;
    std::vector<unsigned int> loc_side_dofs;
    std::vector<unsigned int> loc_edge_dofs;
    unsigned int loc_ele_dof;

    std::shared_ptr<MortarAssemblyBase> mortar_assembly;
};


#endif /* SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_ */
