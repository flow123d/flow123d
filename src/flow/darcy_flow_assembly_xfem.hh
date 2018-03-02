/*
 * darcy_flow_assembly_xfem.hh
 *
 *  Created on: 24, Oct, 2016
 *      Author: pe
 */

#ifndef SRC_FLOW_DARCY_FLOW_ASSEMBLY_XFEM_HH_
#define SRC_FLOW_DARCY_FLOW_ASSEMBLY_XFEM_HH_

#include <memory>
#include "mesh/mesh.h"
#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"
#include "quadrature/quadrature_lib.hh"
#include "quadrature/qmidpoint.hh"
#include "flow/mh_dofhandler.hh"

#include "fem/singularity.hh"
#include "fem/fe_p0_xfem.hh"
#include "fem/fe_rt_xfem.hh"
#include "fem/fe_rt_xfem_single.hh"
#include "quadrature/qxfem.hh"
#include "quadrature/qxfem_factory.hh"

#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
// #include "la/linsys_BDDC.hh"
#include "la/schur.hh"

#include "la/local_to_global_map.hh"
#include "la/local_system.hh"

#include "flow/darcy_flow_mh.hh"
#include "darcy_flow_assembly.hh"

#include "coupling/balance.hh"
#include "flow/mortar_assembly.hh"

template<int dim>
class AssemblyMHXFEM : public AssemblyBase
{
public:
    using AssemblyBase::AssemblyDataPtr;
    
    AssemblyMHXFEM<dim>(AssemblyDataPtr data)
    : quad_(3),
        fe_values_rt_(map_, quad_, fe_rt_,
                update_values | update_gradients | update_JxW_values | update_quadrature_points),

        side_quad_(1),
        fe_p_disc_(0),
        fe_side_values_(map_, side_quad_, fe_p_disc_, update_normal_vectors),

        velocity_interpolation_quad_(0), // veloctiy values in barycenter
        velocity_interpolation_fv_(map_,velocity_interpolation_quad_, fe_rt_, update_values | update_quadrature_points),

        ad_(data),
        sparsity_regular(size(), size()),
        loc_system_(size(), size()),
        loc_system_vb_(2,2)
    {
        unsigned int nsides = RefElement<dim>::n_sides;
        // create regular local sparsity pattern
        sparsity_regular.zeros();
        sparsity_regular.submat(0, 0, nsides, nsides).ones();
        sparsity_regular.diag().ones();
        sparsity_regular.diag(nsides+1).ones();
        sparsity_regular.diag(-nsides-1).ones();
        
        arma::umat sp(2,2);
        // local system 2x2 for vb neighbourings is full matrix
        // this matrix cannot be influenced by any BC (no elimination can take place)
        sp.ones();
        loc_system_vb_.set_sparsity(sp);
    }


    ~AssemblyMHXFEM<dim>() override
    {}

//     LocalSystem& get_local_system() override
//         { return loc_system_;}
    
    void assemble(LocalElementAccessorBase<3> ele_ac) override
    {
        ASSERT_EQ_DBG(ele_ac.dim(), dim);
        
//         DBGVAR(ele_ac.element_accessor().idx());
        setup_local(ele_ac);
        
        if(ele_ac.is_enriched())
            prepare_xfem(ele_ac);
        
        set_dofs_and_bc(ele_ac);
        
        assemble_sides(ele_ac);
        
        assemble_element(ele_ac);
        assemble_source_term(ele_ac);
        
        //mast be last due to overriding xfem fe values
        if(ele_ac.is_enriched())
            assemble_singular_velocity(ele_ac);
        
        assembly_dim_connections(ele_ac);

        if (ad_->balance != nullptr)
            add_fluxes_in_balance_matrix(ele_ac);
        
        if(ele_ac.is_enriched())
        {
            loc_system_.get_matrix().print(cout, "matrix");
            loc_system_.get_rhs().print(cout, "rhs");
        }
        
        ad_->lin_sys->set_local_system(loc_system_);
    }

    void assembly_local_vb(ElementFullIter ele, Neighbour *ngh) override
    {
        ASSERT_LT_DBG(ele->dim(), 3);
        //DebugOut() << "alv " << print_var(this);
        //START_TIMER("Assembly<dim>::assembly_local_vb");
        // compute normal vector to side
        arma::vec3 nv;
        ElementFullIter ele_higher = ad_->mesh->element.full_iter(ngh->side()->element());
        ngh_values_.fe_side_values_.reinit(ele_higher, ngh->side()->el_idx());
        nv = ngh_values_.fe_side_values_.normal_vector(0);

        double value = ad_->sigma.value( ele->centre(), ele->element_accessor()) *
                        2*ad_->conductivity.value( ele->centre(), ele->element_accessor()) *
                        arma::dot(ad_->anisotropy.value( ele->centre(), ele->element_accessor())*nv, nv) *
                        ad_->cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) * // cross-section of higher dim. (2d)
                        ad_->cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) /
                        ad_->cross_section.value( ele->centre(), ele->element_accessor() ) *      // crossection of lower dim.
                        ngh->side()->measure();

        loc_system_vb_.add_value(0,0, -value);
        loc_system_vb_.add_value(0,1,  value);
        loc_system_vb_.add_value(1,0,  value);
        loc_system_vb_.add_value(1,1, -value);
    }


    arma::vec3 make_element_vector_xfem(LocalElementAccessorBase<3> ele_ac){
//         ASSERT_DBG(0).error("Not implemented!");
//         return arma::vec({0,0,0});
        arma::vec3 flux_in_center;
        flux_in_center.zeros();
        ElementFullIter ele = ele_ac.full_iter();
        
        XFEMElementSingularData<dim> * xdata = ele_ac.xfem_data_sing<dim>();
        if(ad_->mh_dh->single_enr) fe_rt_xfem_ = std::make_shared<FE_RT0_XFEM_S<dim,3>>(&fe_rt_,xdata->enrichment_func_vec());
        else fe_rt_xfem_ = std::make_shared<FE_RT0_XFEM<dim,3>>(&fe_rt_,xdata->enrichment_func_vec());
        
        FEValues<dim,3> fv_xfem(map_,velocity_interpolation_quad_, *fe_rt_xfem_, update_values | update_quadrature_points);
        fv_xfem.reinit(ele);
        auto velocity = fv_xfem.vector_view(0);
        
        int dofs[200];
        int ndofs_vel = ele_ac.get_dofs_vel(dofs);
        
        int li = 0;
        for (; li < ndofs_vel; li++) {
            flux_in_center += ad_->mh_dh->mh_solution[dofs[li]] * velocity.value(li,0);
        }
        
        return flux_in_center;
    }
    
    void fix_velocity(LocalElementAccessorBase<3> ele_ac){
        ASSERT_DBG(0).error("Not implemented!");
    }
    
    arma::vec3 make_element_vector(ElementFullIter ele) override
    {
        //START_TIMER("Assembly<dim>::make_element_vector");
//         DBGVAR(ele->index());
        arma::vec3 flux_in_center;
        flux_in_center.zeros();

        //TODO: use LocalElementAccessor
        
        if(ele->xfem_data != nullptr){
            // suppose single proc. so we can create accessor with local ele index
            auto ele_ac = LocalElementAccessorBase<3>(ad_->mh_dh, ele->index());
            flux_in_center = make_element_vector_xfem(ele_ac);
        }
        else{
        
            velocity_interpolation_fv_.reinit(ele);
            for (unsigned int li = 0; li < ele->n_sides(); li++) {
                flux_in_center += ad_->mh_dh->side_flux( *(ele->side( li ) ) )
                            * velocity_interpolation_fv_.vector_view(0).value(li,0);
            }
        }
        
        flux_in_center /= ad_->cross_section.value(ele->centre(), ele->element_accessor() );
        return flux_in_center;
    }

protected:
    static const unsigned int size()
    {
        // sides, 1 for element, edges
        return RefElement<dim>::n_sides + 1 + RefElement<dim>::n_sides;
    }
    
    void prepare_xfem(LocalElementAccessorBase<3> ele_ac){
    
        XFEMElementSingularData<dim> * xdata = ele_ac.xfem_data_sing<dim>();
    
        QXFEMFactory qfact(max_ref_level_[dim]);
        qxfem_ = qfact.create_singular(xdata->sing_vec(), ele_ac.full_iter());
        
    //     qfactory_.gnuplot_refinement<dim>(ele_ac.full_iter(),
    //                                  FilePath("./", FilePath::output_file),
    //                                  *qxfem_,
    //                                  {func});
        
        if(ad_->mh_dh->single_enr) fe_rt_xfem_ = std::make_shared<FE_RT0_XFEM_S<dim,3>>(&fe_rt_,xdata->enrichment_func_vec());
        else fe_rt_xfem_ = std::make_shared<FE_RT0_XFEM<dim,3>>(&fe_rt_,xdata->enrichment_func_vec());
        
        
        fe_values_rt_xfem_ = std::make_shared<FEValues<dim,3>>
                            (map_, *qxfem_, *fe_rt_xfem_, update_values | update_gradients |
                                                        update_JxW_values | update_jacobians |
                                                        update_inverse_jacobians | update_quadrature_points
                                                        | update_divergence);
        
//         fe_p0_xfem_ = std::make_shared<FE_P0_XFEM<dim,3>>(&fe_p_disc_,xdata->enrichment_func_vec());
//         fe_values_p0_xfem_ = std::make_shared<FEValues<dim,3>>
//                             (map_, *qxfem_, *fe_p0_xfem_, update_values |
//                                                         update_JxW_values |
//                                                         update_quadrature_points);
        fe_values_p0_xfem_ = std::make_shared<FEValues<dim,3>>
                            (map_, *qxfem_, fe_p_disc_, update_values |
                                                        update_JxW_values |
                                                        update_quadrature_points);
    }
    
    void setup_local(LocalElementAccessorBase<3> ele_ac){
        
        int dofs[200];
        int ndofs = ele_ac.get_dofs(dofs);
        
        loc_system_.reset(ndofs, ndofs);
        
//         DBGCOUT("####################### DOFS\n");
        for(int i =0; i < ndofs; i++){
//             cout << dofs[i] << " ";
            loc_system_.row_dofs[i] = loc_system_.col_dofs[i] = dofs[i];
        }
//         cout << "\n";
        
        int ndofs_vel = ele_ac.get_dofs_vel(dofs),
            ndofs_pre = ele_ac.get_dofs_press(dofs),
            ndofs_lmb = ele_ac.n_sides();
//         DBGVAR(ndofs_vel);
            
        loc_vel_dofs.resize(ndofs_vel);
        loc_press_dofs.resize(ndofs_pre);
        loc_edge_dofs.resize(ndofs_lmb);
        for(int j =0; j < ndofs_vel; j++)
            loc_vel_dofs[j] = j;
        
        for(int j =0; j < ndofs_pre; j++)
            loc_press_dofs[j] = ndofs_vel + j;
        
        for(int j =0; j < ndofs_lmb; j++)
            loc_edge_dofs[j] = ndofs_vel + ndofs_pre + j;
        
        loc_ele_dof = ndofs_vel;
        
        if(ele_ac.is_enriched())
        {
            int v = ndofs_vel,
                vp = v + ndofs_pre,
                vpl = vp + ndofs_lmb;
            arma::umat sp(ndofs, ndofs);
            sp.zeros();
            sp.diag().ones(); // whole diagonal
            sp.submat(0, 0, v-1, v-1).ones(); // velocity block
            
            sp.submat(v, 0, vp-1, v-1).ones(); // pressure block
            sp.submat(0, v, v-1, vp-1).ones(); // pressure block
            
            sp.submat(vp, 0, vpl-1, v-1).diag().ones(); // lambda block
            sp.submat(0, vp, v-1, vpl-1).diag().ones(); // lambda block
            
            if(ndofs > vpl){
                sp.submat(vpl, 0, ndofs-1, v-1).ones(); // lambda_w block
                sp.submat(0, vpl, v-1, ndofs-1).ones(); // lambda_w block
            }
//             sp.print("SP_enriched");
            loc_system_.set_sparsity(sp);
        }
        else
            loc_system_.set_sparsity(sparsity_regular);
    }

    void set_dofs_and_bc(LocalElementAccessorBase<3> ele_ac){
        
        ASSERT_DBG(ele_ac.dim() == dim);
        
        //set global dof for element (pressure)
//         loc_system_.row_dofs[loc_ele_dof] = loc_system_.col_dofs[loc_ele_dof] = ele_ac.ele_row();
        
        //shortcuts
        const unsigned int nsides = ele_ac.n_sides();
        LinSys *ls = ad_->lin_sys;
        
        Boundary *bcd;
        unsigned int side_row, edge_row;
        
        dirichlet_edge.resize(nsides);
        for (unsigned int i = 0; i < nsides; i++) {

            side_row = loc_vel_dofs[i];    //local
            edge_row = loc_edge_dofs[i];    //local
            
//             auto c = ele_ac.side(i)->centre();
//             if(ele_ac.dim() == 2)
//             DBGCOUT(<< "ele " << ele_ac.ele_global_idx() << " s " << loc_system_.row_dofs[side_row] << " e " << loc_system_.row_dofs[edge_row]
//                     << " c " << c(0) << " " << c(1) << " " << c(2) << "\n");
            
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
            
//                     DBGCOUT(<< "[" << loc_system_.row_dofs[edge_row] << ", " << loc_system_.row_dofs[edge_row]
//                             << "] mat: " << -bcd->element()->measure() * bc_sigma * cross_section
//                             << " rhs: " << (bc_flux - bc_sigma * bc_pressure) * bcd->element()->measure() * cross_section
//                             << "\n");
                    dirichlet_edge[i] = 2;  // to be skipped in LMH source assembly
                    loc_system_.add_value(edge_row, edge_row,
                                            -bcd->element()->measure() * bc_sigma * cross_section,
                                            (bc_flux - bc_sigma * bc_pressure) * bcd->element()->measure() * cross_section);
                }
                else if (type==DarcyMH::EqData::seepage) {
                    ad_->is_linear=false;

                    unsigned int loc_edge_idx = bcd->bc_ele_idx_;
                    char & switch_dirichlet = ad_->bc_switch_dirichlet[loc_edge_idx];
                    double bc_pressure = ad_->bc_switch_pressure.value(b_ele.centre(), b_ele);
                    double bc_flux = -ad_->bc_flux.value(b_ele.centre(), b_ele);
                    double side_flux = bc_flux * bcd->element()->measure() * cross_section;

                    // ** Update BC type. **
                    if (switch_dirichlet) {
                        // check and possibly switch to flux BC
                        // The switch raise error on the corresponding edge row.
                        // Magnitude of the error is abs(solution_flux - side_flux).
                        ASSERT_DBG(ad_->mh_dh->rows_ds->is_local(ele_ac.side_row(i)))(ele_ac.side_row(i));
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
                        ASSERT_DBG(ad_->mh_dh->rows_ds->is_local(ele_ac.edge_row(i)))(ele_ac.edge_row(i));
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
                            loc_system_.set_solution(edge_row,bc_pressure, -1);
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
                    ASSERT_DBG(ad_->mh_dh->rows_ds->is_local(ele_ac.edge_row(i)))(ele_ac.edge_row(i));
                    unsigned int loc_edge_row = ele_ac.edge_local_row(i);
                    double & solution_head = ls->get_solution_array()[loc_edge_row];

                    // Force Robin type during the first iteration of the unsteady case.
                    if (solution_head > bc_switch_pressure  || ad_->force_bc_switch) {
                        // Robin BC
                        //DebugOut().fmt("x: {}, robin, bcp: {}\n", b_ele.centre()[0], bc_pressure);
                        loc_system_.add_value(edge_row, edge_row,
                                                -bcd->element()->measure() * bc_sigma * cross_section,
                                                bcd->element()->measure() * cross_section * (bc_flux - bc_sigma * bc_pressure)  );
                    } else {
                        // Neumann BC
                        //DebugOut().fmt("x: {}, neuman, q: {}  bcq: {}\n", b_ele.centre()[0], bc_switch_pressure, bc_pressure);
                        double bc_total_flux = bc_flux + bc_sigma*(bc_switch_pressure - bc_pressure);
                        
                        loc_system_.add_value(edge_row, bc_total_flux * bcd->element()->measure() * cross_section);
                    }
                } 
                else {
                    xprintf(UsrErr, "BC type not supported.\n");
                }
            }
            loc_system_.add_value(side_row, edge_row, 1.0);
            loc_system_.add_value(edge_row, side_row, 1.0);
            
            
//             if(ad_->mh_dh->single_enr)
//             if(ad_->mh_dh->enrich_velocity && ele_ac.is_enriched()){
//                 assemble_enriched_side_edge(ele_ac, i);
//             }
        }
    }
    
    
    void assemble_enriched_side_edge(LocalElementAccessorBase<3> ele_ac, unsigned int local_side){
    }
    
    
    void assemble_sides(LocalElementAccessorBase<3> ele_ac) override
    {
        double cs = ad_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor());
        double conduct =  ad_->conductivity.value(ele_ac.centre(), ele_ac.element_accessor());
        double scale = 1 / cs /conduct;
        
        if(ele_ac.is_enriched())
            assemble_sides_scale(ele_ac, scale, *fe_values_rt_xfem_);
        else
            assemble_sides_scale(ele_ac, scale, fe_values_rt_);
    }
    
    void assemble_sides_scale(LocalElementAccessorBase<3> ele_ac, double scale, FEValues<dim,3> & fe_values)
    {
        arma::vec3 &gravity_vec = ad_->gravity_vec_;
        
        ElementFullIter ele =ele_ac.full_iter();
        fe_values.reinit(ele);
        unsigned int ndofs = loc_vel_dofs.size();
        unsigned int qsize = fe_values.get_quadrature()->size();

        auto velocity = fe_values.vector_view(0);
        
        for (unsigned int k=0; k<qsize; k++)
            for (unsigned int i=0; i<ndofs; i++){
                double rhs_val =
                        arma::dot(gravity_vec,velocity.value(i,k))
                        * fe_values.JxW(k);
                loc_system_.add_value(i, rhs_val);
                
                for (unsigned int j=0; j<ndofs; j++){
                    double mat_val = 
                        arma::dot(velocity.value(i,k), //TODO: compute anisotropy before
                                    (ad_->anisotropy.value(ele_ac.centre(), ele_ac.element_accessor() )).i()
                                        * velocity.value(j,k))
                        * scale * fe_values.JxW(k);
                    
                    loc_system_.add_value(i, j, mat_val);
                }
            }
        
//         // assemble matrix for weights in BDDCML
//         // approximation to diagonal of 
//         // S = -C - B*inv(A)*B'
//         // as 
//         // diag(S) ~ - diag(C) - 1./diag(A)
//         // the weights form a partition of unity to average a discontinuous solution from neighbouring subdomains
//         // to a continuous one
//         // it is important to scale the effect - if conductivity is low for one subdomain and high for the other,
//         // trust more the one with low conductivity - it will be closer to the truth than an arithmetic average
//         if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
//             const arma::mat& local_matrix = loc_system_.get_matrix();
//             for(unsigned int i=0; i < ndofs; i++) {
//                 double val_side =  local_matrix(i,i);
//                 double val_edge =  -1./local_matrix(i,i);
// 
//                 unsigned int side_row = loc_system_.row_dofs[loc_side_dofs[i]];
//                 unsigned int edge_row = loc_system_.row_dofs[loc_edge_dofs[i]];
//                 static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( side_row, val_side );
//                 static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( edge_row, val_edge );
//             }
//         }
    }
    
    
    void assemble_element(LocalElementAccessorBase<3> ele_ac){        
        
        if(ele_ac.is_enriched()){
                assemble_element(ele_ac, *fe_values_rt_xfem_, *fe_values_p0_xfem_);
        }
        else 
        {
            for(unsigned int side = 0; side < ele_ac.n_sides(); side++){
                loc_system_.add_value(loc_ele_dof, loc_vel_dofs[side], -1.0);
                loc_system_.add_value(loc_vel_dofs[side], loc_ele_dof, -1.0);
            }
        }
        
//         if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
//             double val_ele =  1.;
//             static_cast<LinSys_BDDC*>(ad_->lin_sys)->
//                             diagonal_weights_set_value( loc_system_.row_dofs[loc_ele_dof], val_ele );
//         }
    }
    
    void assemble_element(LocalElementAccessorBase<3> ele_ac,
                          FEValues<dim,3>& fv_vel, FEValues<dim,3>& fv_press){
        
        ElementFullIter ele = ele_ac.full_iter();
        fv_vel.reinit(ele);
        fv_press.reinit(ele);
        
        unsigned int ndofs_vel = loc_vel_dofs.size();
        unsigned int ndofs_press = loc_press_dofs.size();
        unsigned int qsize = qxfem_->size();

        for (unsigned int k=0; k<qsize; k++)
            for (unsigned int i=0; i<ndofs_vel; i++){
//             for (unsigned int i=0; i<ele_ac.n_sides(); i++){
                for (unsigned int j=0; j<ndofs_press; j++){
                    double mat_val = 
                        - fv_press.shape_value(j,k)
                        * fv_vel.shape_divergence(i,k)
                        * fv_vel.JxW(k);
                    loc_system_.add_value(loc_vel_dofs[i],loc_press_dofs[j] , mat_val, 0.0);
                    loc_system_.add_value(loc_press_dofs[j], loc_vel_dofs[i] , mat_val, 0.0);
                }
            }
        
//         loc_system_.add_value(ele_ac.n_sides(), loc_press_dofs[0], 2*M_PI, 0.0);
//         loc_system_.add_value(loc_press_dofs[0], ele_ac.n_sides(), 2*M_PI, 0.0);
        
//         for (unsigned int i=0; i<ndofs_vel; i++)
//             for (unsigned int j=0; j<ndofs_press; j++)
//                 DBGCOUT(<< i << " " << j << " " << loc_system_.get_matrix()(loc_vel_dofs[i],loc_press_dofs[j]) 
//                         << " " << ele_ac.full_iter()->measure()  << "  " << fv_vel.determinant(0) << "\n");
    }
    
    void assembly_dim_connections(LocalElementAccessorBase<3> ele_ac){
        //D, E',E block: compatible connections: element-edge
        ElementFullIter ele = ele_ac.full_iter();
        
        // no Neighbours => nothing to asssemble here
        if(ele->n_neighs_vb == 0) return;
        
        int ele_row = ele_ac.ele_row();
        Neighbour *ngh;

        //DebugOut() << "adc " << print_var(this) << print_var(side_quad_.size());
        for (unsigned int i = 0; i < ele->n_neighs_vb; i++) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            ngh = ele_ac.full_iter()->neigh_vb[i];
            loc_system_vb_.reset();
            loc_system_vb_.row_dofs[0] = loc_system_vb_.col_dofs[0] = ele_row;
            loc_system_vb_.row_dofs[1] = loc_system_vb_.col_dofs[1] = ad_->mh_dh->row_4_edge[ ngh->edge_idx() ];

            assembly_local_vb(ele, ngh);

            ad_->lin_sys->set_local_system(loc_system_vb_);

            // update matrix for weights in BDDCML
            if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
               int ind = loc_system_vb_.row_dofs[1];
               // there is -value on diagonal in block C!
               double new_val = loc_system_vb_.get_matrix()(0,0);
               static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( ind, new_val );
            }
        }
    }
    
    void add_fluxes_in_balance_matrix(LocalElementAccessorBase<3> ele_ac){

        for (unsigned int i = 0; i < ele_ac.n_sides(); i++) {
            Boundary* bcd = ele_ac.side(i)->cond();

            if (bcd) {
                /*
                    DebugOut().fmt("add_flux: {} {} {} {}\n",
                            ad_->mh_dh->el_ds->myp(),
                            ele_ac.ele_global_idx(),
                            ad_->local_boundary_index,
                            ele_ac.side_row(i));
                 */
                ad_->balance->add_flux_matrix_values(ad_->water_balance_idx, ad_->local_boundary_index,
                                                     {(int)(ele_ac.side_row(i))}, {1});
                ++(ad_->local_boundary_index);
            }
        }
    }
    
    void assemble_singular_velocity(LocalElementAccessorBase<3> ele_ac){
        XFEMElementSingularData<dim> * xd = ele_ac.xfem_data_sing<dim>();
        ElementFullIter ele = ele_ac.full_iter();
        
        //as long as pressure is not enriched and is P0
        ASSERT_DBG(! ad_->mh_dh->enrich_pressure);

        double cs = ad_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor());
        double val;
        
        for(unsigned int w=0; w < xd->n_enrichments(); w++){
            if(xd->enrichment_intersects(w)){
                auto sing = static_pointer_cast<Singularity<dim-2>>(xd->enrichment_func(w));

                auto quad = xd->sing_quadrature(w);
                fv_rt_sing_ = std::make_shared<FEValues<dim,3>>
                                (map_, quad, *fe_rt_xfem_, update_values);

                fv_rt_sing_->reinit(ele);
                auto velocity = fv_rt_sing_->vector_view(0);

                unsigned int loc_sing_dof = loc_edge_dofs[0] + loc_edge_dofs.size() + w;
    //             DBGVAR(loc_sing_dof);
                
                // robin like condition with sigma
                // well lagrange multiplier test function is 1
                
                // local part of the effective_surface in the element
                double effective_surface = 0;
                for(unsigned int q=0; q < quad.size();q++)
                    effective_surface += quad.weight(q);
                
                effective_surface = cs * effective_surface;
                
                for(unsigned int q=0; q < quad.size();q++){
                    // outer normal is opposite to distance vector
                    arma::vec n = - sing->geometry().dist_vector(quad.real_point(q));
                    n = n / arma::norm(n,2);
                    
                    // assembly well boundary integral
                    
                    for (int i=0; i < loc_vel_dofs.size(); i++){
                        val = arma::dot(velocity.value(i,q),n)
                            * quad.weight(q);
                        
                        loc_system_.add_value(loc_vel_dofs[i], loc_sing_dof, val, 0.0);
                        loc_system_.add_value(loc_sing_dof, loc_vel_dofs[i], val, 0.0);
                    }
                }
                
//                 loc_system_.add_value(loc_sing_dof, loc_sing_dof,
//                                     - effective_surface * sing->sigma(),
//                                     - effective_surface * sing->sigma() * sing->pressure());
                
                // Bw integral p1-lambda_w            
                double sing_lagrange_val = 1;
//                 if(quad.size() == sing->q_points().size()) val = sing->circumference() * sing->sigma();
//                 else{
//                     val = 0;
//                     for(unsigned int q=0; q < quad.size();q++)
//                         val += quad.weight(q);
//                     val = sing_lagrange_val * sing->sigma() * val;
//                 }
                val = sing_lagrange_val * sing->sigma() * effective_surface;
                
                DBGVAR(val);
                int sing_row = ele_ac.sing_row(w);
                int ele1d_row = ad_->mh_dh->row_4_el[xd->intersection_ele_global_idx()];
                DBGVAR(sing_row);
                DBGVAR(ele1d_row);
                ad_->lin_sys->mat_set_value(sing_row, ele1d_row, val);
                ad_->lin_sys->mat_set_value(ele1d_row, sing_row, val);
                ad_->lin_sys->mat_set_value(sing_row, sing_row, -val);
                ad_->lin_sys->mat_set_value(ele1d_row, ele1d_row, -val);
            }
        }
    }
    
//     void assemble_singular_velocity(LocalElementAccessorBase<3> ele_ac){
//         XFEMElementSingularData<dim> * xd = ele_ac.xfem_data_sing<dim>();
//         ElementFullIter ele = ele_ac.full_iter();
//         
//         //as long as pressure is not enriched and is P0
//         ASSERT_DBG(! ad_->mh_dh->enrich_pressure);
// 
//         double cs = ad_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor());
//         double val;
//         
//         for(unsigned int w=0; w < xd->n_enrichments(); w++){
//             if(xd->enrichment_intersects(w)){
//                 auto sing = static_pointer_cast<Singularity<dim-2>>(xd->enrichment_func(w));
// 
//                 auto quad = xd->sing_quadrature(w);
//                 fv_rt_sing_ = std::make_shared<FEValues<dim,3>>
//                                 (map_, quad, *fe_rt_xfem_, update_values);
// 
//                 fv_rt_sing_->reinit(ele);
//                 auto velocity = fv_rt_sing_->vector_view(0);
// 
//                 unsigned int loc_sing_dof = loc_edge_dofs[0] + loc_edge_dofs.size() + w;
//     //             DBGVAR(loc_sing_dof);
//                 
//                 // robin like condition with sigma
//                 // well lagrange multiplier test function is 1
//                 
//                 // local part of the effective_surface in the element
//                 double effective_surface = 0;
//                 for(unsigned int q=0; q < quad.size();q++)
//                     effective_surface += quad.weight(q);
//                 
//                 effective_surface = cs * effective_surface;
//                 
//                 for(unsigned int q=0; q < quad.size();q++){
//                     // outer normal is opposite to distance vector
//                     arma::vec n = - sing->geometry().dist_vector(quad.real_point(q));
//                     n = n / arma::norm(n,2);
//                     
//                     // assembly well boundary integral
//                     
//                     for (int i=0; i < loc_vel_dofs.size(); i++){
//                         val = arma::dot(velocity.value(i,q),n)
//                             * quad.weight(q);
//                         
//                         loc_system_.add_value(loc_vel_dofs[i], loc_sing_dof, val, 0.0);
//                         loc_system_.add_value(loc_sing_dof, loc_vel_dofs[i], val, 0.0);
//                     }
//                 }
//                 
//                 loc_system_.add_value(loc_sing_dof, loc_sing_dof,
//                                     - effective_surface * sing->sigma(),
//                                     - effective_surface * sing->sigma() * sing->pressure());
//             }
//         }
//     }
    
    // assembly volume integrals
    FE_RT0<dim,3> fe_rt_;
    MappingP1<dim,3> map_;
    QGauss<dim> quad_;
    FEValues<dim,3> fe_values_rt_;
    
    NeighSideValues<dim<3?dim:2> ngh_values_;

    // assembly face integrals (BC)
    QGauss<dim-1> side_quad_;
    FE_P_disc<dim,3> fe_p_disc_;
    FESideValues<dim,3> fe_side_values_;

    // Interpolation of velocity into barycenters
    QGauss<dim> velocity_interpolation_quad_;
    FEValues<dim,3> velocity_interpolation_fv_;

    // data shared by assemblers of different dimension
    AssemblyDataPtr ad_;
    std::vector<unsigned int> dirichlet_edge;

    arma::umat sparsity_regular;
    LocalSystem loc_system_;
    LocalSystem loc_system_vb_;
    std::vector<unsigned int> loc_edge_dofs;
    unsigned int loc_ele_dof;
    
    
    int dof_tmp[100];
    
    // assembly volume integrals
    const std::vector<unsigned int> max_ref_level_ = {1, 12, 12, 5};
    shared_ptr<QXFEM<dim,3>> qxfem_;
    
    shared_ptr<FiniteElementEnriched<dim,3>> fe_rt_xfem_;
    shared_ptr<FEValues<dim,3>> fe_values_rt_xfem_;
    
    shared_ptr<FiniteElementEnriched<dim,3>> fe_p0_xfem_;
    shared_ptr<FEValues<dim,3>> fe_values_p0_xfem_;
    
    shared_ptr<FEValues<dim,3>> fv_rt_sing_;
    
    std::vector<unsigned int> loc_vel_dofs;
    std::vector<unsigned int> loc_press_dofs;
};

template<> inline void AssemblyMHXFEM<1>::prepare_xfem(LocalElementAccessorBase<3> ele_ac){}
template<> inline void AssemblyMHXFEM<1>::assemble_singular_velocity(LocalElementAccessorBase<3> ele_ac){}


// template<> inline
// void AssemblyMHXFEM<2>::assemble_singular_velocity(LocalElementAccessorBase<3> ele_ac){
// 
//     ElementFullIter ele = ele_ac.full_iter();
//     
//     XFEMElementSingularData * xd = ele_ac.xfem_data_sing();
//     double sing_lagrange_val;
//     int sing_row;
//     int nvals = loc_vel_dofs.size();
//     double val[nvals];
//     int vel_dofs[nvals];
//     
//     for(unsigned int w=0; w < xd->n_enrichments(); w++){
//         if(xd->is_singularity_inside(w)){
//             sing_row = ele_ac.sing_row(w);
//             auto quad = xd->sing_quadrature(w);
//             fv_rt_sing_ = std::make_shared<FEValues<2,3>> 
//                             (map_, quad, *fe_rt_xfem_, update_values);
// 
//             fv_rt_sing_->reinit(ele);
//             auto sing = static_pointer_cast<Singularity0D<3>>(xd->enrichment_func(w));
//             sing_lagrange_val = 1;
//             for (int i=0; i < nvals; i++){
//                 val[i] = 0;
//                 vel_dofs[i] = loc_system_.row_dofs[i];
// //                 double sum = 0;
//                 for(unsigned int q=0; q < quad.size();q++){
//                     arma::vec n = sing->center() - quad.real_point(q);
//                     n = n / arma::norm(n,2);
//                     val[i] += sing_lagrange_val
//                                 * arma::dot(fv_rt_sing_->shape_vector(i,q),n)
//                                 * quad.weight(q);
// //                     sum += quad.weight(q);
//                 }
// //                 DBGVAR(sum);
// //                 DBGVAR(sing->circumference());
//             }
//             
//             ad_->lin_sys->mat_set_values(1, &sing_row, nvals, vel_dofs, val);
//             ad_->lin_sys->mat_set_values(nvals, vel_dofs, 1, &sing_row, val);
//         }
//     }
// }

// template<> inline
// void AssemblyMHXFEM<2>::assemble_singular_velocity(LocalElementAccessorBase<3> ele_ac){
// 
//     ElementFullIter ele = ele_ac.full_iter();
//     XFEMElementSingularData * xd = ele_ac.xfem_data_sing();
//     int sing_row;
//     double sing_lagrange_val = 1.0;
//     
//     //as long as pressure is not enriched and is P0
//     ASSERT_DBG(! ad_->mh_dh->enrich_pressure);
//     int ele1d_row = ad_->mh_dh->row_4_el[xd->intersection_ele_global_idx()];
//     
//     int nvals = loc_vel_dofs.size();
//     double val;
//     double vals[nvals];
//     int vel_dofs[nvals];
//     
//     for(unsigned int w=0; w < xd->n_enrichments(); w++){
//         if(xd->is_singularity_inside(w)){
//             auto quad = xd->sing_quadrature(w);
//             fv_rt_sing_ = std::make_shared<FEValues<2,3>> 
//                             (map_, quad, *fe_rt_xfem_, update_values);
// 
//             fv_rt_sing_->reinit(ele);
//             auto sing = static_pointer_cast<Singularity0D<3>>(xd->enrichment_func(w));
//             sing_row = ele_ac.sing_row(w);
//             
//             // Bw integral p2-(v.n)
//             for (int i=0; i < nvals; i++){
//                 vals[i] = 0;
//                 vel_dofs[i] = loc_system_.row_dofs[i];
// //                 double sum = 0;
//                 for(unsigned int q=0; q < quad.size();q++){
//                     arma::vec n = sing->center() - quad.real_point(q);
//                     n = n / arma::norm(n,2);
//                     vals[i] += sing_lagrange_val
//                                 * arma::dot(fv_rt_sing_->shape_vector(i,q),n)
//                                 * quad.weight(q);
// //                     sum += quad.weight(q);
//                 }
// //                 DBGVAR(sum);
// //                 DBGVAR(sing->circumference());
// //                 DBGVAR(val[i]);
// //                 loc_system_.add_value(loc_ele_dof, loc_vel_dofs[i], val, 0.0);
// //                 loc_system_.add_value(loc_vel_dofs[i], loc_ele_dof, val, 0.0);
//             }
//             ad_->lin_sys->mat_set_values(1, &sing_row, nvals, vel_dofs, vals);
//             ad_->lin_sys->mat_set_values(nvals, vel_dofs, 1, &sing_row, vals);
//             
//             
//             // Bw integral p1-lambda_w            
//             if(quad.size() == sing->q_points().size()) val = sing->circumference() * sing->sigma();
//             else{
//                 val = 0;
//                 for(unsigned int q=0; q < quad.size();q++)
//                     val += quad.weight(q);
//                 val = sing_lagrange_val * sing->sigma() * val;
//             }
//             
//             DBGVAR(val);
//             DBGVAR(sing_row);
//             DBGVAR(ele1d_row);
//             ad_->lin_sys->mat_set_value(sing_row, ele1d_row, val);
//             ad_->lin_sys->mat_set_value(ele1d_row, sing_row, val);
//             ad_->lin_sys->mat_set_value(sing_row, sing_row, -val);
//             ad_->lin_sys->mat_set_value(ele1d_row, ele1d_row, -val);
//         }
//     }
// }

// template<> inline
// void AssemblyMHXFEM<2>::assemble_singular_velocity(LocalElementAccessorBase<3> ele_ac){
// 
// //     ElementFullIter ele = ele_ac.full_iter();
//     XFEMElementSingularData<2> * xd = ele_ac.xfem_data_sing<2>();
//     
//     //as long as pressure is not enriched and is P0
//     ASSERT_DBG(! ad_->mh_dh->enrich_pressure);
// //     double temp = 1.0;
// //     int ele1d_row = ad_->mh_dh->row_4_el[xd->intersection_ele_global_idx()];
//     
// //     int nvals = loc_vel_dofs.size();
// //     double val;
//     
//     
//     for(unsigned int w=0; w < xd->n_enrichments(); w++){
//         if(xd->enrichment_intersects(w)){
// //             auto quad = xd->sing_quadrature(w);
// //             fv_rt_sing_ = std::make_shared<FEValues<2,3>> 
// //                             (map_, quad, *fe_rt_xfem_, update_values);
// 
// //             fv_rt_sing_->reinit(ele);
//             auto sing = static_pointer_cast<Singularity0D>(xd->enrichment_func(w));
// //             temp = 1.0 / sing->sigma();
//             
// //             vector<double> sum(nvals,0);
// //             arma::mat matt(nvals, nvals+1);
//             
// //             for(unsigned int q=0; q < quad.size();q++){
// //                 arma::vec n = sing->center() - quad.real_point(q);
// //                 n = n / arma::norm(n,2);
// //                 for (int i=0; i < nvals; i++){
//                     // int B_w 1/sigma * (u.n)*(v.n)
// //                     for (int j=0; j < nvals; j++){
// //                         val = temp 
// //                             * arma::dot(fv_rt_sing_->shape_vector(i,q),n)
// //                             * arma::dot(fv_rt_sing_->shape_vector(j,q),n)
// //                             * quad.weight(q);
// //                         loc_system_.add_value(loc_vel_dofs[i], loc_vel_dofs[j], val, 0.0);
// // //                         loc_system_.add_value(loc_vel_dofs[j], loc_vel_dofs[i], val, 0.0);
// //                         matt(i,j) += val;
// //                         
// //                     }
//                     //dirichlet on the RHS
// //                     val = - sing->pressure() * arma::dot(fv_rt_sing_->shape_vector(i,q),n) * quad.weight(q);
// //                     if(i==3) DBGVAR(val);
// //                     sum[i] += val;
// //                     loc_system_.add_value(loc_vel_dofs[i], loc_vel_dofs[i], 0.0, val);
// //                     matt(i,nvals) = val;
// //                     sum[i] += val;
// //                 }
// //             }
// //             DBGCOUT(<< "\n" << matt);
// //             DBGVAR(sum[0]);
// //             DBGVAR(sum[1]);
// //             DBGVAR(sum[2]);
// //             DBGVAR(sum[3]);
// 
//             unsigned int loc_sing_dof = loc_edge_dofs[0] + loc_edge_dofs.size() + w;
//             unsigned int loc_enr_vel_dof = loc_vel_dofs[fe_rt_xfem_->n_regular_dofs()] + w;
//             DBGVAR(loc_sing_dof);
//             DBGVAR(loc_enr_vel_dof);
//             // robin like condition with sigma
//             // well lagrange multiplier test function is 1
//             double effective_surface = sing->geometry().effective_surface();
//             double sing_flow = arma::norm(sing->vector(xd->sing_quadrature(w).real_point(0)),2) * effective_surface;
//             DBGVAR(sing_flow);
//             loc_system_.add_value(loc_enr_vel_dof, loc_sing_dof, sing_flow, 0);
//             loc_system_.add_value(loc_sing_dof, loc_enr_vel_dof, sing_flow, 0);
//             
//             loc_system_.add_value(loc_sing_dof, loc_sing_dof,
//                                   - effective_surface * sing->sigma(),
//                                   - effective_surface * sing->sigma() * sing->pressure());
//         }
//     }
// }

// template<> inline
// void AssemblyMHXFEM<2>::assemble_enriched_side_edge(LocalElementAccessorBase<3> ele_ac, unsigned int local_side){
//         ElementFullIter ele = ele_ac.full_iter();
// //         ele->node[0]->point().print(cout, "point 0");
// //         ele->node[1]->point().print(cout, "point 1");
// //         ele->node[2]->point().print(cout, "point 2");
//         DBGVAR(local_side);
//         
//         //Simply create normal vector.
//         QGauss<1> auxq(1);
//         auto fv_side = std::make_shared<FESideValues<2,3>>(map_, auxq, *fe_rt_xfem_, update_normal_vectors);
//         fv_side->reinit(ele, local_side);
//         
// //         fv_side->normal_vector(0).print(cout,"normal");
//         
//         const uint qsize=100; // 1d quadrature on side
//         QXFEM<2,3> qside_xfem(QMidpoint(qsize), local_side, *ele->permutation_idx_); // mapped side quadrature to 2d coords
//         for(unsigned int q=0; q < qsize; q++){   // map to real coords
//             arma::vec real_point = map_.project_unit_to_real(RefElement<2>::local_to_bary(qside_xfem.point(q)),map_.element_map(*ele));
//             qside_xfem.set_real_point(q,real_point);
//         }
//         
//         auto fv_xfem = std::make_shared<FEValues<2,3>>(map_, qside_xfem, *fe_rt_xfem_, update_values);
//         fv_xfem->reinit(ele);
//         
//         for(unsigned int j=fe_rt_xfem_->n_regular_dofs(); j<fe_rt_xfem_->n_dofs(); j++){
// //             side_row = loc_system_.row_dofs[loc_vel_dofs[j]];
// //             DBGVAR(j);
// //             DBGVAR(loc_vel_dofs[j]);
// //             DBGVAR(loc_edge_dofs[local_side]);
//                     
//             double sum_val = 0;
//             double side_measure = ele->side(local_side)->measure();
//             for(unsigned int q=0; q < qsize; q++){
// //                 auto qp = qside_xfem.real_point(q);
// //                 cout << qp(0) << " " << qp(1) << " " << qp(2) << "\n";
// //                 auto fv = fv_xfem->shape_vector(j,q);
// //                 cout << fv(0) << " " << fv(1) << " " << fv(2) << "\n";
//                 double val = arma::dot(fv_xfem->shape_vector(j,q),fv_side->normal_vector(0))
//                       // this makes JxW on the triangle side:
//                       * qside_xfem.weight(q)
//                       * side_measure;
//                       
// //                 ad_->lin_sys->mat_set_value(side_row, edge_row, val);
// //                 ad_->lin_sys->mat_set_value(edge_row, side_row, val);
//                    loc_system_.add_value(loc_vel_dofs[j], loc_edge_dofs[local_side], val, 0.0);
//                    loc_system_.add_value(loc_edge_dofs[local_side], loc_vel_dofs[j], val, 0.0);
//                 sum_val += val;
//             }
//             DBGVAR(sum_val);
// //             DBGVAR(side_measure);
//         }        
//     }

template<> inline
void AssemblyMHXFEM<2>::assemble_enriched_side_edge(LocalElementAccessorBase<3> ele_ac, unsigned int local_side){
        ElementFullIter ele = ele_ac.full_iter();
        DBGVAR(local_side);
        
        XFEMElementSingularData<2> * xdata = ele_ac.xfem_data_sing<2>();
        
        for(unsigned int j=fe_rt_xfem_->n_regular_dofs(); j<fe_rt_xfem_->n_dofs(); j++){
//             side_row = loc_system_.row_dofs[loc_vel_dofs[j]];
//             DBGVAR(j);
//             DBGVAR(loc_vel_dofs[j]);
//             DBGVAR(loc_edge_dofs[local_side]);
            //Simply create normal vector.
            QGauss<1> auxq(1);
            auto fv_side = std::make_shared<FESideValues<2,3>>(map_, auxq, *fe_rt_xfem_, update_normal_vectors);
            fv_side->reinit(ele, local_side);
            
            QXFEMFactory qfact(max_ref_level_[1]);
            auto qside_xfem = qfact.create_side_singular(xdata->sing_vec(),
                                                             ele_ac.full_iter(), local_side);
            auto fv_xfem = std::make_shared<FEValues<2,3>>(map_, *qside_xfem, *fe_rt_xfem_, update_values);
            fv_xfem->reinit(ele);
            auto velocity = fv_xfem->vector_view(0);
            
            double sum_val = 0;
//             double side_measure = ele->side(local_side)->measure();
            for(unsigned int q=0; q < qside_xfem->size(); q++){
//                 auto qp = qside_xfem.real_point(q);
//                 cout << qp(0) << " " << qp(1) << " " << qp(2) << "\n";
//                 auto fv = velocity.value(j,q);
//                 cout << fv(0) << " " << fv(1) << " " << fv(2) << "\n";
                double val = arma::dot(velocity.value(j,q),fv_side->normal_vector(0))
                           * qside_xfem->JxW(q);
                      // this makes JxW on the triangle side:
//                       * qside_xfem->weight(q)
//                       * side_measure;
                      
//                 ad_->lin_sys->mat_set_value(side_row, edge_row, val);
//                 ad_->lin_sys->mat_set_value(edge_row, side_row, val);
                   loc_system_.add_value(loc_vel_dofs[j], loc_edge_dofs[local_side], val, 0.0);
                   loc_system_.add_value(loc_edge_dofs[local_side], loc_vel_dofs[j], val, 0.0);
                sum_val += val;
            }
            DBGVAR(sum_val);
//             DBGVAR(side_measure);
        }        
    }
    
template<> inline
void AssemblyMHXFEM<3>::assemble_enriched_side_edge(LocalElementAccessorBase<3> ele_ac, unsigned int local_side){
        ElementFullIter ele = ele_ac.full_iter();
        DBGVAR(local_side);
//         const double side_measure = ele->side(local_side)->measure();
//         DBGVAR(side_measure);
        
        XFEMElementSingularData<3> * xdata = ele_ac.xfem_data_sing<3>();

        //Simply create normal vector.
        QGauss<2> auxq(1);
        auto fv_side = std::make_shared<FESideValues<3,3>>(map_, auxq, *fe_rt_xfem_, update_normal_vectors);
        fv_side->reinit(ele, local_side);
        
        QXFEMFactory qfact(max_ref_level_[2]);
        auto qside_xfem = qfact.create_side_singular(xdata->sing_vec(),
                                                            ele_ac.full_iter(), local_side);
        auto fv_xfem = std::make_shared<FEValues<3,3>>(map_, *qside_xfem, *fe_rt_xfem_, update_values);
        fv_xfem->reinit(ele);
        auto velocity = fv_xfem->vector_view(0);
            
        for(unsigned int j=fe_rt_xfem_->n_regular_dofs(); j<fe_rt_xfem_->n_dofs(); j++){
//             DBGVAR(j);
//             DBGVAR(loc_vel_dofs[j]);
//             DBGVAR(loc_edge_dofs[local_side]);

            double val, sum_val = 0;
            
            for(unsigned int q=0; q < qside_xfem->size(); q++){
//                 auto qp = qside_xfem.real_point(q);
//                 cout << qp(0) << " " << qp(1) << " " << qp(2) << "\n";
//                 auto fv = velocity.value(j,q);
//                 cout << fv(0) << " " << fv(1) << " " << fv(2) << "\n";
                val = arma::dot(velocity.value(j,q),fv_side->normal_vector(0))
                    * qside_xfem->JxW(q);   
                    // this makes JxW on the triangle side:
//                     * qside_xfem->weight(q)
//                     * side_measure;
                    
                loc_system_.add_value(loc_vel_dofs[j], loc_edge_dofs[local_side], val, 0.0);
                loc_system_.add_value(loc_edge_dofs[local_side], loc_vel_dofs[j], val, 0.0);
                sum_val += val;
            }
            DBGVAR(sum_val);
        }        
    }
    
#endif /* SRC_FLOW_DARCY_FLOW_ASSEMBLY_XFEM_HH_ */
