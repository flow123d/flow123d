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
        fe_side_values_(map_, side_quad_, fe_p_disc_, update_normal_vectors),

        velocity_interpolation_quad_(0), // veloctiy values in barycenter
        velocity_interpolation_fv_(map_,velocity_interpolation_quad_, fe_rt_, update_values | update_quadrature_points),

        ad_(data),
        loc_system_(size(), size())
    {
    }


    ~AssemblyMHXFEM<dim>() override
    {}

    LocalSystem& get_local_system() override
        { return loc_system_;}
    
    void assemble(LocalElementAccessorBase<3> ele_ac) override
    {
        ASSERT_DBG(ele_ac.dim() == dim);
        
//         unsigned int ndofs = ele_ac.n_dofs();
        setup_local(ele_ac);
        
        if(ele_ac.is_enriched() && !ele_ac.xfem_data_pointer()->is_complement())
            prepare_xfem(ele_ac);
        
        set_dofs_and_bc(ele_ac);
        
        assemble_sides(ele_ac);
        
        assemble_element(ele_ac);
        assemble_source_term(ele_ac);
        
        //mast be last due to overriding xfem fe values
        if(ele_ac.is_enriched() && !ele_ac.xfem_data_pointer()->is_complement())
            assemble_singular_velocity(ele_ac);
        
        if(ele_ac.is_enriched())
            loc_system_.get_matrix().print(cout, "matrix");
        
        loc_system_.fix_diagonal();
    }

    void assembly_local_vb(double *local_vb,  ElementFullIter ele, Neighbour *ngh) override
    {
        //START_TIMER("Assembly<dim>::assembly_local_vb");
        // compute normal vector to side
        arma::vec3 nv;
        ElementFullIter ele_higher = ad_->mesh->element.full_iter(ngh->side()->element());
        fe_side_values_.reinit(ele_higher, ngh->side()->el_idx());
        nv = fe_side_values_.normal_vector(0);

        double value = ad_->sigma.value( ele->centre(), ele->element_accessor()) *
                        2*ad_->conductivity.value( ele->centre(), ele->element_accessor()) *
                        arma::dot(ad_->anisotropy.value( ele->centre(), ele->element_accessor())*nv, nv) *
                        ad_->cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) * // cross-section of higher dim. (2d)
                        ad_->cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) /
                        ad_->cross_section.value( ele->centre(), ele->element_accessor() ) *      // crossection of lower dim.
                        ngh->side()->measure();

        local_vb[0] = -value;   local_vb[1] = value;
        local_vb[2] = value;    local_vb[3] = -value;
    }


    arma::vec3 make_element_vector_xfem(LocalElementAccessorBase<3> ele_ac){
        ASSERT_DBG(0).error("Not implemented!");
        return arma::vec({0,0,0});
    }
    
    arma::vec3 make_element_vector(ElementFullIter ele) override
    {
        //START_TIMER("Assembly<dim>::make_element_vector");
        DBGVAR(ele->index());
        arma::vec3 flux_in_center;
        flux_in_center.zeros();

        //TODO: use LocalElementAccessor
        
        if(ele->xfem_data != nullptr && ! ele->xfem_data->is_complement()){
            // suppose single proc. so we can create accessor with local ele index
            auto ele_ac = LocalElementAccessorBase<3>(ad_->mh_dh, ele->index());
            flux_in_center = make_element_vector_xfem(ele_ac);
        }
        else{
        
            velocity_interpolation_fv_.reinit(ele);
            for (unsigned int li = 0; li < ele->n_sides(); li++) {
                flux_in_center += ad_->mh_dh->side_flux( *(ele->side( li ) ) )
                            * velocity_interpolation_fv_.shape_vector(li,0);
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
    
    void prepare_xfem(LocalElementAccessorBase<3> ele_ac)
    { }
    
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
            ndofs_pre = ele_ac.get_dofs_press(dofs);
//         DBGVAR(ndofs_vel);
            
        loc_vel_dofs.resize(ndofs_vel);
        loc_press_dofs.resize(ndofs_pre);
        loc_edge_dofs.resize(ele_ac.n_sides());
        for(int j =0; j < ndofs_vel; j++)
            loc_vel_dofs[j] = j;
        
        for(int j =0; j < ndofs_pre; j++)
            loc_press_dofs[j] = ndofs_vel + j;
        
        for(int j =0; j < ele_ac.n_sides(); j++)
            loc_edge_dofs[j] = ndofs_vel + ndofs_pre + j;
        
        loc_ele_dof = ndofs_vel;
        
//         //velocity
//         int ndofs = ele_ac.get_dofs_vel(dof_tmp);
//         for (unsigned int i = 0; i < ndofs; i++) {
//             loc_system_.col_dofs[loc_vel_dofs[i]] = dof_tmp[i];
//         }
//         
//         //pressure
//         ndofs = ele_ac.get_dofs_press(dof_tmp);
//         for (unsigned int i = 0; i < ndofs; i++) {
//             loc_system_.col_dofs[loc_press_dofs[i]] = dof_tmp[i];
//         }
            
    }

    void set_dofs_and_bc(LocalElementAccessorBase<3> ele_ac){
        
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
                    loc_system_.set_solution(edge_row,bc_pressure);
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
                    loc_system_.set_value(edge_row, edge_row,
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
                            loc_system_.set_solution(edge_row,bc_pressure);
                            dirichlet_edge[i] = 1;
                        } else {
                            //DebugOut()("x: {}, neuman, q: {}  bcq: {}\n", b_ele.centre()[0], side_flux, bc_flux);
                            loc_system_.set_value(edge_row, side_row, 1.0, side_flux);
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
                        loc_system_.set_value(edge_row, edge_row,
                                                -bcd->element()->measure() * bc_sigma * cross_section,
                                                bcd->element()->measure() * cross_section * (bc_flux - bc_sigma * bc_pressure)  );
                    } else {
                        // Neumann BC
                        //DebugOut().fmt("x: {}, neuman, q: {}  bcq: {}\n", b_ele.centre()[0], bc_switch_pressure, bc_pressure);
                        double bc_total_flux = bc_flux + bc_sigma*(bc_switch_pressure - bc_pressure);
                        
                        loc_system_.set_value(edge_row, side_row,
                                                1.0,
                                                bc_total_flux * bcd->element()->measure() * cross_section);
                    }
                } 
                else {
                    xprintf(UsrErr, "BC type not supported.\n");
                }
            }
            loc_system_.set_mat_values({side_row}, {edge_row}, {1.0});
            loc_system_.set_mat_values({edge_row}, {side_row}, {1.0});
            
            
            if(ad_->mh_dh->enrich_velocity && ele_ac.is_enriched() && ! ele_ac.xfem_data_pointer()->is_complement()){
                assemble_enriched_side_edge(ele_ac, i);
            }
        }
    }
    
    
    void assemble_enriched_side_edge(LocalElementAccessorBase<3> ele_ac, unsigned int local_side){
    }
    
    
    void assemble_sides(LocalElementAccessorBase<3> ele_ac) override
    {
        double cs = ad_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor());
        double conduct =  ad_->conductivity.value(ele_ac.centre(), ele_ac.element_accessor());
        double scale = 1 / cs /conduct;
        
        if(ele_ac.is_enriched() && !ele_ac.xfem_data_pointer()->is_complement())
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

        for (unsigned int k=0; k<qsize; k++)
            for (unsigned int i=0; i<ndofs; i++){
                double rhs_val =
                        arma::dot(gravity_vec,fe_values.shape_vector(i,k))
                        * fe_values.JxW(k);
                loc_system_.add_value(i,i , 0.0, rhs_val);
                
                for (unsigned int j=0; j<ndofs; j++){
                    double mat_val = 
                        arma::dot(fe_values.shape_vector(i,k), //TODO: compute anisotropy before
                                    (ad_->anisotropy.value(ele_ac.centre(), ele_ac.element_accessor() )).i()
                                        * fe_values.shape_vector(j,k))
                        * scale * fe_values.JxW(k);
                    
                    loc_system_.add_value(i,j , mat_val, 0.0);
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
        
        if(ele_ac.is_enriched() && !ele_ac.xfem_data_pointer()->is_complement()){
                assemble_element(ele_ac, *fe_values_rt_xfem_, *fe_values_p0_xfem_);
        }
        else 
        {
            for(unsigned int side = 0; side < ele_ac.n_sides(); side++){
                loc_system_.set_mat_values({loc_ele_dof}, {loc_vel_dofs[side]}, {-1.0});
                loc_system_.set_mat_values({loc_vel_dofs[side]}, {loc_ele_dof}, {-1.0});
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
                for (unsigned int j=0; j<ndofs_press; j++){
                    double mat_val = 
                        - fv_press.shape_value(j,k)
                        * fv_vel.shape_divergence(i,k)
                        * fv_vel.JxW(k);
                    loc_system_.add_value(loc_vel_dofs[i],loc_press_dofs[j] , mat_val, 0.0);
                    loc_system_.add_value(loc_press_dofs[j], loc_vel_dofs[i] , mat_val, 0.0);
                }
            }
        
//         for (unsigned int i=0; i<ndofs_vel; i++)
//             for (unsigned int j=0; j<ndofs_press; j++)
//                 DBGCOUT(<< i << " " << j << " " << loc_system_.get_matrix()(loc_vel_dofs[i],loc_press_dofs[j]) 
//                         << " " << ele_ac.full_iter()->measure()  << "  " << fv_vel.determinant(0) << "\n");
    }
    
    void assemble_singular_velocity(LocalElementAccessorBase<3> ele_ac){
    }
    
    // assembly volume integrals
    FE_RT0<dim,3> fe_rt_;
    MappingP1<dim,3> map_;
    QGauss<dim> quad_;
    FEValues<dim,3> fe_values_rt_;

    // assembly face integrals (BC)
    QGauss<dim-1> side_quad_;
    FE_P_disc<0,dim,3> fe_p_disc_;
    FESideValues<dim,3> fe_side_values_;

    // Interpolation of velocity into barycenters
    QGauss<dim> velocity_interpolation_quad_;
    FEValues<dim,3> velocity_interpolation_fv_;

    // data shared by assemblers of different dimension
    AssemblyDataPtr ad_;
    std::vector<unsigned int> dirichlet_edge;

    LocalSystem loc_system_;
    std::vector<unsigned int> loc_edge_dofs;
    unsigned int loc_ele_dof;
    
    
    int dof_tmp[100];
    
    // assembly volume integrals
    QXFEMFactory<dim,3> qfactory_;
    shared_ptr<QXFEM<dim,3>> qxfem_;
    
    shared_ptr<FiniteElementEnriched<dim,3>> fe_rt_xfem_;
    shared_ptr<FEValues<dim,3>> fe_values_rt_xfem_;
    shared_ptr<FESideValues<dim,3>> fv_side_xfem_;
    
    shared_ptr<FiniteElementEnriched<dim,3>> fe_p0_xfem_;
    shared_ptr<FEValues<dim,3>> fe_values_p0_xfem_;
    
    shared_ptr<FEValues<dim,3>> fv_rt_sing_;
    
    std::vector<unsigned int> loc_vel_dofs;
    std::vector<unsigned int> loc_press_dofs;
};














template<> inline
void AssemblyMHXFEM<2>::prepare_xfem(LocalElementAccessorBase<3> ele_ac){
        
    XFEMElementSingularData * xdata = ele_ac.xfem_data_sing();
    
    std::shared_ptr<Singularity0D<3>> func = std::static_pointer_cast<Singularity0D<3>>(xdata->enrichment_func(0));
    
    qxfem_ = qfactory_.create_singular({func}, ele_ac.full_iter());
//     qfactory_.gnuplot_refinement(ele_ac.full_iter(),
//                                  FilePath("./", FilePath::output_file),
//                                  *qxfem_,
//                                  {func});
    
    if(ad_->mh_dh->single_enr) fe_rt_xfem_ = std::make_shared<FE_RT0_XFEM_S<2,3>>(&fe_rt_,xdata->enrichment_func_vec());
    else fe_rt_xfem_ = std::make_shared<FE_RT0_XFEM<2,3>>(&fe_rt_,xdata->enrichment_func_vec());
    
    
    fe_values_rt_xfem_ = std::make_shared<FEValues<2,3>> 
                        (map_, *qxfem_, *fe_rt_xfem_, update_values | update_gradients |
                                                    update_JxW_values | update_jacobians |
                                                    update_inverse_jacobians | update_quadrature_points 
                                                    | update_divergence);
    
    fe_p0_xfem_ = std::make_shared<FE_P0_XFEM<2,3>>(&fe_p_disc_,xdata->enrichment_func_vec());
    fe_values_p0_xfem_ = std::make_shared<FEValues<2,3>> 
                        (map_, *qxfem_, *fe_p0_xfem_, update_values |
                                                    update_JxW_values |
                                                    update_quadrature_points);
}


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

template<> inline
void AssemblyMHXFEM<2>::assemble_singular_velocity(LocalElementAccessorBase<3> ele_ac){

    ElementFullIter ele = ele_ac.full_iter();
    XFEMElementSingularData * xd = ele_ac.xfem_data_sing();
    
    //as long as pressure is not enriched and is P0
    ASSERT_DBG(! ad_->mh_dh->enrich_pressure);
    double temp = 1.0;
    int ele1d_row = ad_->mh_dh->row_4_el[xd->intersection_ele_global_idx()];
    
    int nvals = loc_vel_dofs.size();
    double val;
    
    for(unsigned int w=0; w < xd->n_enrichments(); w++){
        if(xd->is_singularity_inside(w)){
            auto quad = xd->sing_quadrature(w);
            fv_rt_sing_ = std::make_shared<FEValues<2,3>> 
                            (map_, quad, *fe_rt_xfem_, update_values);

            fv_rt_sing_->reinit(ele);
            auto sing = static_pointer_cast<Singularity0D<3>>(xd->enrichment_func(w));
            temp = 1.0 / sing->sigma();
            
//             vector<double> sum(nvals,0);
            for(unsigned int q=0; q < quad.size();q++){
                arma::vec n = sing->center() - quad.real_point(q);
                n = n / arma::norm(n,2);
                for (int i=0; i < nvals; i++){
                    // int B_w 1/sigma * (u.n)*(v.n)
                    for (int j=0; j < nvals; j++){
                        val = temp 
                            * arma::dot(fv_rt_sing_->shape_vector(i,q),n)
                            * arma::dot(fv_rt_sing_->shape_vector(j,q),n)
                            * quad.weight(q);
                        loc_system_.add_value(loc_vel_dofs[i], loc_vel_dofs[j], val, 0.0);
                        loc_system_.add_value(loc_vel_dofs[j], loc_vel_dofs[i], val, 0.0);
                    }
                    
//                     DBGVAR(val);
                    val = - sing->pressure() * arma::dot(fv_rt_sing_->shape_vector(i,q),n) * quad.weight(q);
                    loc_system_.add_value(loc_vel_dofs[i], loc_vel_dofs[i], 0.0, val);
                    
//                     sum[i] += val;
//                     ad_->lin_sys->mat_set_value(loc_system_.row_dofs[loc_vel_dofs[i]], ele1d_row, val);
//                     ad_->lin_sys->mat_set_value(ele1d_row, loc_system_.row_dofs[loc_vel_dofs[i]], val);
                }
            }
//             DBGVAR(sum[0]);
//             DBGVAR(sum[1]);
//             DBGVAR(sum[2]);
        }
    }
}


template<> inline
void AssemblyMHXFEM<2>::assemble_enriched_side_edge(LocalElementAccessorBase<3> ele_ac, unsigned int local_side){
        ElementFullIter ele = ele_ac.full_iter();
        double val;
        DBGVAR(local_side);
        
//         int side_row, edge_row;
//         edge_row = loc_system_.row_dofs[loc_edge_dofs[local_side]];
        
//         QGauss<1> qside(29);
        QGauss<1> qside(1);
        const uint qsize=100;
        qside.resize(qsize);
        double qweight = 1.0/qsize;
        for(unsigned int q=0; q<qsize; q++){
            qside.set_point(q, arma::vec({0.5*qweight + q*qweight}));
            qside.set_weight(q,qweight);
        }
        
        
        
        auto fv_side = std::make_shared<FESideValues<2,3>>(map_, qside, *fe_rt_xfem_,
                                                            update_normal_vectors
                                                            | update_quadrature_points);
        fv_side->reinit(ele, local_side);
        
        QXFEM<2,3> qside_xfem;
        qside_xfem.resize(qside.size());
        for(unsigned int q=0; q < qside.size(); q++){
            //TODO: use mapping for quad from side to element
            arma::vec unit_p = ele->project_point(fv_side->point(q));
            qside_xfem.set_point(q, RefElement<2>::bary_to_local(unit_p));
            qside_xfem.set_real_point(q,fv_side->point(q));
//             fv_side->point(q).print(cout, "side point");
            qside_xfem.set_weight(q,qside.weight(q));
        }
        
        auto fv_xfem = std::make_shared<FEValues<2,3>>(map_, qside_xfem, *fe_rt_xfem_,
                                                        update_values);
        fv_xfem->reinit(ele);
        
        for(unsigned int j=fe_rt_xfem_->n_regular_dofs(); j<fe_rt_xfem_->n_dofs(); j++){
//             side_row = loc_system_.row_dofs[loc_vel_dofs[j]];
//             DBGVAR(j);
//             DBGVAR(loc_vel_dofs[j]);
//             DBGVAR(loc_edge_dofs[local_side]);
                    
            double sum_val = 0;
            double side_measure = ele->side(local_side)->measure();
            for(unsigned int q=0; q < qside.size(); q++){
                val = arma::dot(fv_xfem->shape_vector(j,q),fv_side->normal_vector(q))
                      // this makes JxW on the triangle side:
                      * qside_xfem.weight(q)
                      * side_measure;
                      
//                 ad_->lin_sys->mat_set_value(side_row, edge_row, val);
//                 ad_->lin_sys->mat_set_value(edge_row, side_row, val);
                   loc_system_.add_value(loc_vel_dofs[j], loc_edge_dofs[local_side], val, 0.0);
                   loc_system_.add_value(loc_edge_dofs[local_side], loc_vel_dofs[j], val, 0.0);
                sum_val += val;
            }
            DBGVAR(sum_val);
        }        
        

        
//         fv_side_xfem_ = std::make_shared<FESideValues<2,3>>(map_, side_quad_, *fe_rt_xfem_,
//                                                             update_values
//                                                             | update_normal_vectors);
//         fv_side_xfem_->reinit(ele, local_side);
// 
//         for(unsigned int j=fe_rt_xfem_->n_regular_dofs(); j<fe_rt_xfem_->n_dofs(); j++){
//              side_row = loc_system_.row_dofs[loc_vel_dofs[j]];
//             //suppose one point quadrature at the moment
//             val = arma::dot(fv_side_xfem_->shape_vector(j,0),fv_side_xfem_->normal_vector(0))
//                     * ele->side(local_side)->measure();
//             ad_->lin_sys->mat_set_value(side_row, edge_row, val);
//             ad_->lin_sys->mat_set_value(edge_row, side_row, val);
//         }
        
    }
    

template<> inline
arma::vec3 AssemblyMHXFEM<2>::make_element_vector_xfem(LocalElementAccessorBase<3> ele_ac){
    arma::vec3 flux_in_center;
    flux_in_center.zeros();
    ElementFullIter ele = ele_ac.full_iter();
    
    XFEMElementSingularData * xdata = ele_ac.xfem_data_sing();
    fe_rt_xfem_ = std::make_shared<FE_RT0_XFEM<2,3>>(&fe_rt_,xdata->enrichment_func_vec());
    
    FEValues<2,3> fv_xfem(map_,velocity_interpolation_quad_, *fe_rt_xfem_, update_values | update_quadrature_points);
    fv_xfem.reinit(ele);
    
    unsigned int li = 0;
    for (; li < ele->n_sides(); li++) {
        flux_in_center += ad_->mh_dh->side_flux( *(ele->side( li ) ) )
                        * fv_xfem.shape_vector(li,0);
    }
    
    for (unsigned int w = 0; w < xdata->n_enrichments(); w++){
        std::vector<int>& wdofs = xdata->global_enriched_dofs()[Quantity::velocity][w];
        for (unsigned int j = 0; j < wdofs.size(); j++, li++) {
            flux_in_center += ad_->mh_dh->mh_solution[wdofs[li]]
                            * fv_xfem.shape_vector(li,0);
        }
    }
    
    return flux_in_center;
}    
    
#endif /* SRC_FLOW_DARCY_FLOW_ASSEMBLY_XFEM_HH_ */