/*
 * mortar_assembly.cc
 *
 *  Created on: Feb 22, 2017
 *      Author: jb
 */

#include "flow/mortar_assembly.hh"
#include "quadrature/quadrature_lib.hh"
#include "quadrature/intersection_quadrature.hh"
#include "la/linsys.hh"
//#include "mesh/intersection.hh"
#include "intersection/mixed_mesh_intersections.hh"
#define ARMA_USE_CXX11
#include <armadillo>

P0_CouplingAssembler::P0_CouplingAssembler(AssemblyDataPtr data)
: MortarAssemblyBase(data),
  tensor_average_(16),
  col_average_(4),
  quadrature_(*(data->mesh)),
  slave_ac_(data->mh_dh)
{
    isec_data_list.reserve(30);


    for(uint row_dim=0; row_dim<4; row_dim ++)
        for(uint col_dim=0; col_dim<4; col_dim ++) {
            arma::vec row_avg = arma::vec(row_dim+1);
            row_avg.fill(1.0 / (row_dim+1));
            arma::vec col_avg = arma::vec(col_dim+1);
            col_avg.fill(1.0 / (col_dim+1));
            arma::mat avg = row_avg * col_avg.t(); // tensor product
            tensor_average(row_dim, col_dim) = avg;
            col_average_[row_dim] = row_avg;
        }
}


void P0_CouplingAssembler::pressure_diff(LocalElementAccessorBase<3> ele_ac, double delta) {
    isec_data_list.push_back(IsecData());
    IsecData &i_data = isec_data_list.back();

    i_data.dim= ele_ac.dim();
    i_data.delta = delta;
    i_data.dofs.resize(ele_ac.n_sides());
    i_data.vel_dofs.resize(ele_ac.n_sides());
    i_data.dirichlet_dofs.resize(ele_ac.n_sides());
    i_data.dirichlet_sol.resize(ele_ac.n_sides());
    i_data.n_dirichlet=0;
    i_data.ele_z_coord_=ele_ac.centre()[2];

    for(unsigned int i_side=0; i_side < ele_ac.n_sides(); i_side++ ) {
        i_data.dofs[i_side]=ele_ac.edge_row(i_side);
        i_data.vel_dofs[i_side] = ele_ac.side_row(i_side);
        //i_data.z_sides[i_side]=ele_ac.side(i_side)->centre()[2];
        //DebugOut().fmt("edge: {} {}", i_side, ele_ac.edge_row(i_side));
        Boundary * bcd = ele_ac.full_iter()->side(i_side)->cond();
        if (bcd) {
            ElementAccessor<3> b_ele = bcd->element_accessor();
            auto type = (DarcyMH::EqData::BC_Type)data_->bc_type.value(b_ele.centre(), b_ele);
            //DebugOut().fmt("bcd id: {} sidx: {} type: {}\n", ele->id(), i_side, type);
            if (type == DarcyMH::EqData::dirichlet) {
                //DebugOut().fmt("Dirichlet: {}\n", ele->index());
                double bc_pressure = data_->bc_pressure.value(b_ele.centre(), b_ele);
                i_data.dirichlet_dofs[i_data.n_dirichlet] = i_side;
                i_data.dirichlet_sol[i_data.n_dirichlet] = bc_pressure;
                i_data.n_dirichlet++;
            }
        }
    }

}





/**
 * Works well but there is large error next to the boundary.
 */
void P0_CouplingAssembler::assembly(LocalElementAccessorBase<3> master_ac)
{


    // on the intersection element we consider
    // * intersection dofs for master and slave
    //   those are dofs of the space into which we interpolate
    //   base functions from individual master and slave elements
    //   For the master dofs both are usualy eqivalent.
    // * original dofs - for master same as intersection dofs, for slave
    //   all dofs of slave elements

    // form list of intersection dofs, in this case pressures in barycenters
    // but we do not use those form MH system in order to allow second schur somplement. We rather map these
    // dofs to pressure traces, i.e. we use averages of traces as barycentric values.

    /**
     *  TODO:
     *  - pass through the master and all slaves and collect global dofs , bcd, solution.
     *    I.e. call Nx pressure_diff not NxNx.
     *
     *  - Is it safe to have duplicate rows in local_system?
     *  - Is it better to have more smaller local system then single big one?
     *  - use one big or more smaller local systems to set.
     */


    if (master_ac.dim() > 2) return; // supported only for 1D and 2D master elements
    auto &isec_list = mixed_mesh_.element_intersections_[master_ac.ele_global_idx()];
    if (isec_list.size() == 0) return; // skip empty masters

    //slave_ac_.setup(master_ac);

    ElementFullIter ele = master_ac.full_iter();
    arma::vec3 ele_centre = ele->centre();
    double m_sigma = data_->sigma.value( ele_centre, ele->element_accessor());
    double m_conductivity = data_->conductivity.value( ele_centre, ele->element_accessor());
    double m_crossection = data_->cross_section.value( ele_centre, ele->element_accessor() );

    double master_sigma = 2*m_sigma*m_conductivity *
                    2/ m_crossection;
                    /**
                     * ?? How to deal with anisotropy ??
                     * 3d-2d : compute nv of 2d triangle
                     * 2d-2d : interpret as 2d-1d-2d, should be symmetric master-slave
                     * 2d-1d : nv is tangent to 2d and normal to 1d
                    arma::dot(data_->anisotropy.value( ele_centre, ele->element_accessor())*nv, nv)
                    */







    isec_data_list.clear();
    double cs_sqr_avg = 0.0;
    double isec_sum = 0.0;
    uint i = 0;
    for(; i < isec_list.size(); ++i) {
        bool non_zero = quadrature_.reinit(isec_list[i].second);
        slave_ac_.reinit( quadrature_.slave_idx() );
        if (slave_ac_.dim() == master_ac.dim()) break;
        if (! non_zero) continue; // skip quadratures close to zero

        double cs = data_->cross_section.value(slave_ac_.full_iter()->centre(), slave_ac_.full_iter()->element_accessor());
        double isec_measure = quadrature_.measure();
        //DebugOut() << print_var(cs) << print_var(isec_measure);
        cs_sqr_avg += cs*cs*isec_measure;
        isec_sum += isec_measure;
        //DebugOut().fmt("Assembly23: {} {} {} ", ele->id(), slave_ac_.full_iter()->id(), isec_measure);
        pressure_diff(slave_ac_, isec_measure);
    }
    element_isec = isec_sum;
    if ( ! (slave_ac_.dim() == 2 && master_ac.dim() ==2 ) ) {
        if( fabs(isec_sum - ele->measure()) > 1E-5) {
            string out;
            for(auto & isec : isec_list) {
                slave_ac_.reinit(isec.second->bulk_ele_idx());
                out += fmt::format(" {}", slave_ac_.full_iter()->id());
            }

            double diff = (isec_sum - ele->measure())/ele->measure();
            WarningOut() << print_var(diff)
                    << print_var(ele->id())
                    << endl
                    << out;

        }
    }
    pressure_diff(master_ac, -isec_sum);

    //DebugOut().fmt( "cs2: {} d0: {}", cs_sqr_avg, delta_0);
    master_sigma = master_sigma * (cs_sqr_avg / isec_sum)
            / isec_sum;


    add_to_linsys(master_sigma);


    // 2d-2d
    //DebugOut() << "2d-2d";
    if (i < isec_list.size()) {
        isec_data_list.clear();
        isec_sum = 0.0;
        for(; i < isec_list.size(); ++i) {
                quadrature_.reinit(isec_list[i].second);
                slave_ac_.reinit( quadrature_.slave_idx() );
                double isec_measure = quadrature_.measure();
                isec_sum += isec_measure;
                //DebugOut().fmt("Assembly22: {} {} {}", ele->id(), slave_ac_.full_iter()->id(), isec_measure);
                pressure_diff(slave_ac_, isec_measure);
        }
        pressure_diff(master_ac, -isec_sum);

        master_sigma = 100 * m_conductivity/ m_crossection / isec_sum;

        add_to_linsys(master_sigma);
    }
}

void P0_CouplingAssembler::add_to_linsys(double scale)
{
    // rows
     for(IsecData &row_ele : isec_data_list) {
         //columns
         for(IsecData &col_ele : isec_data_list) {


             double s =  -scale * row_ele.delta * col_ele.delta;
             //DebugOut().fmt("s: {} {} {} {}", s, scale, row_ele.delta, col_ele.delta);
             product_ = s * tensor_average(row_ele.dim, col_ele.dim);
             //product_ = s * row_ele.values * col_ele.values;

             loc_system_.reset(row_ele.dofs, col_ele.dofs);

             for(uint i=0; i< row_ele.n_dirichlet; i++)
                 loc_system_.set_solution_row(row_ele.dirichlet_dofs[i], row_ele.dirichlet_sol[i], -1.0);
             for(uint i=0; i< col_ele.n_dirichlet; i++) loc_system_.set_solution_col(col_ele.dirichlet_dofs[i], col_ele.dirichlet_sol[i]);
             //ASSERT( arma::norm(product_,2) == 0.0 );
             loc_system_.set_matrix(product_);
             // Must have z-coords for every side, can not use averaging vector
             loc_system_.set_rhs( -s * col_average_[row_ele.dim] * col_ele.ele_z_coord_ ); // TODO: should be row_ele.ele_z_coord_, however it possibly doesnt matter due to symmetry
             loc_system_.eliminate_solution();

             if (fix_velocity_flag) {
                 this->fix_velocity_local(row_ele, col_ele);
             } else {
                 data_->lin_sys->set_local_system(loc_system_);
             }
         }
     }
}


void P0_CouplingAssembler::fix_velocity_local(const IsecData &row_ele, const IsecData & col_ele)
{

    uint n_rows = row_ele.vel_dofs.n_rows;
    uint n_cols = col_ele.dofs.n_rows;
    arma::vec pressure(n_cols);
    arma::vec add_velocity(n_rows);
    double * solution = data_->lin_sys->get_solution_array();
    for(uint icol=0; icol < n_cols; icol++ ) pressure[icol] = solution[col_ele.dofs[icol]];
    add_velocity =  loc_system_.get_matrix() * pressure - loc_system_.get_rhs();
    //DebugOut() << "fix_velocity\n" << pressure << add_velocity;
    for(uint irow=0; irow < n_rows; irow++ ) solution[row_ele.vel_dofs[irow]] += add_velocity[irow] ;
}












P1_CouplingAssembler::P1_CouplingAssembler(AssemblyDataPtr data)
: MortarAssemblyBase(data),
  slave_ac_(data->mh_dh),
  master_ac_(data->mh_dh)
{
    /// TODO: Use RefElement numberings to make it independent.
    P1_for_P1e[0] = arma::mat({{1}});
    P1_for_P1e[1] = arma::eye(2,2);

    //{{1, 1, -1}, {1, -1, 1}, {-1, 1, 1}}
    P1_for_P1e[2] = arma::ones(3,3);
    P1_for_P1e[2].at(0,2)=-1;
    P1_for_P1e[2].at(1,1)=-1;
    P1_for_P1e[2].at(2,0)=-1;

    //{{1, 1, 1, -1}, {1, 1, -1, 1}, {1, -1, 1, 1}, {-1, 1, 1, 1}}
    P1_for_P1e[3] = arma::ones(4,4);
    P1_for_P1e[3].at(0,3)=-2;
    P1_for_P1e[3].at(1,2)=-2;
    P1_for_P1e[3].at(2,1)=-2;
    P1_for_P1e[3].at(3,0)=-2;

}


void P1_CouplingAssembler::add_sides(LocalElementAccessorBase<3> ele_ac, unsigned int shift, arma::uvec &dofs)
{
    for(unsigned int i_side=0; i_side < ele_ac.n_sides(); i_side++ ) {
        uint loc_dof = shift+i_side;
        dofs[loc_dof] =  ele_ac.edge_row(i_side);
        Boundary * bcd = ele_ac.full_iter()->side(i_side)->cond();

        if (bcd) {
            ElementAccessor<3> b_ele = bcd->element_accessor();
            auto type = (DarcyMH::EqData::BC_Type)data_->bc_type.value(b_ele.centre(), b_ele);
            //DebugOut().fmt("bcd id: {} sidx: {} type: {}\n", ele->id(), i_side, type);
            if (type == DarcyMH::EqData::dirichlet) {
                double bc_pressure = data_->bc_pressure.value(b_ele.centre(), b_ele);
                loc_system_.set_solution(loc_dof, bc_pressure, -1.0);
            }
        }
    }
}


template<uint qdim, uint master_dim, uint slave_dim>
void P1_CouplingAssembler::isec_assembly(
        double master_sigma,
        const IntersectionLocal<master_dim, slave_dim> &il,
        std::array<uint, qdim+1 > subdiv)
{
    // mappings from quadrature ref. el (bary) to intersection ref. elements (local coords)
    arma::mat master_map(master_dim, qdim + 1);
    arma::mat slave_map(slave_dim, qdim+1);
    for(uint i_col=0; i_col < qdim+1; i_col++) {
        uint ip = subdiv[i_col];
        master_map.col(i_col) = il[ip].comp_coords();
        slave_map.col(i_col) = il[ip].bulk_coords();
    }

    uint n_dofs=master_ac_.n_sides() + slave_ac_.n_sides();
    arma::uvec dofs(n_dofs);
    loc_system_.reset(n_dofs, n_dofs);
    add_sides(master_ac_, 0, dofs);
    add_sides(slave_ac_, master_ac_.n_sides(), dofs);
    loc_system_.set_dofs(dofs, dofs);

    QGauss<qdim> q(2);
    arma::vec values(n_dofs);

    element_isec = 0.0;
    for(uint ip=0; ip < q.size(); ip++) {
        double il_measure = il.compute_measure();
        arma::vec bary_qp =  RefElement<qdim>::local_to_bary(q.point(ip));

        arma::vec qp_master = RefElement<master_dim>::local_to_bary(master_map * bary_qp);
        arma::vec qp_slave = RefElement<slave_dim>::local_to_bary(slave_map * bary_qp);
        double JxW = q.weight(ip) * il_measure * master_jac;
        uint shift = 0;
        values.rows(shift, shift + master_dim) = P1_for_P1e[master_dim] * qp_master;
        shift=master_dim + 1;
        values.rows(shift, shift + slave_dim ) = -P1_for_P1e[slave_dim] * qp_slave;

        double scale = -master_sigma * JxW;
        arma::mat x = scale * ( values * values.t() );

        //DebugOut() << print_var(scale);
        //DebugOut() << print_var(values);
        //DebugOut() << x;
        loc_system_.get_matrix() += x;

        element_isec+= il_measure;
    }

    loc_system_.eliminate_solution();
    if (fix_velocity_flag) {
        //this->fix_velocity_local(row_ele, col_ele);
    } else {
        data_->lin_sys->set_local_system(loc_system_);
    }


}



/**
 * P1 connection of different dimensions
 *
 * - 20.11. 2014 - very poor convergence, big error in pressure even at internal part of the fracture
 */

void P1_CouplingAssembler::assembly(LocalElementAccessorBase<3> master_ac) {

    if (master_ac.dim() > 2) return; // supported only for 1D and 2D master elements
    auto &isec_list = mixed_mesh_.element_intersections_[master_ac.ele_global_idx()];
    if (isec_list.size() == 0) return; // skip empty masters

    //slave_ac_.setup(master_ac);

    ElementFullIter master_ele = master_ac.full_iter();
    master_ac_.reinit(master_ac.ele_global_idx());
    arma::vec3 ele_centre = master_ele->centre();
    double m_sigma = data_->sigma.value( ele_centre, master_ele->element_accessor());
    double m_conductivity = data_->conductivity.value( ele_centre, master_ele->element_accessor());
    double m_crossection = data_->cross_section.value( ele_centre, master_ele->element_accessor() );

    double master_sigma = 2*m_sigma*m_conductivity *
                    2/ m_crossection;
                    /**
                     * ?? How to deal with anisotropy ??
                     * 3d-2d : compute nv of 2d triangle
                     * 2d-2d : interpret as 2d-1d-2d, should be symmetric master-slave
                     * 2d-1d : nv is tangent to 2d and normal to 1d
                    arma::dot(data_->anisotropy.value( ele_centre, ele->element_accessor())*nv, nv)
                    */

    master_jac = master_ele->measure() * master_ele->dim();

    for(auto &isec_pair : isec_list) {
        IntersectionLocalBase *isec = isec_pair.second;
        ASSERT_EQ_DBG(master_ac.ele_global_idx(), isec->component_ele_idx());

        //ElementFullIter slave_ele = data_->mesh->element(isec->bulk_ele_idx());
        slave_ac_.reinit(isec->bulk_ele_idx());


        if (typeid(*isec) == typeid(IntersectionLocal<1,2>)) {
            auto il = static_cast<const IntersectionLocal<1,2> *>(isec);
            ASSERT_EQ_DBG( il->size(), 2);
            DebugOut() << "il12";

            this->isec_assembly<1,1,2>(master_sigma, *il, {0,1});

        } else
        if (typeid(*isec) == typeid(IntersectionLocal<2,2>)) {
            auto il = static_cast<const IntersectionLocal<2,2> *>(isec);
            ASSERT_EQ_DBG( il->size(), 2);
            DebugOut() << "il22";

            this->isec_assembly<1,2,2>(master_sigma, *il, {0,1});

        } else if (typeid(*isec) == typeid(IntersectionLocal<2,3>)) {
            auto il = static_cast<const IntersectionLocal<2,3> *>(isec);
            if (il->size() <= 2) continue; // skip degenerated intersections
            const uint master_dim = 2;
            const uint slave_dim = 3;
            const uint qdim = 2;

            uint n_dofs=master_ac_.n_sides() + slave_ac_.n_sides();
            arma::uvec dofs(n_dofs);
            loc_system_.reset(n_dofs, n_dofs);
            add_sides(master_ac_, 0, dofs);
            add_sides(slave_ac_, master_ac_.n_sides(), dofs);
            loc_system_.set_dofs(dofs, dofs);

            double isec_area = 0;
            double isec_measure = il->compute_measure();

            //DebugOut() << "subdivision, n: " << il->size();

            for(uint i_vtx=2; i_vtx< il->size();  i_vtx++) {
                //this->isec_assembly<2,2,3>(master_sigma, *il, {0, i_vtx1, i_vtx1+1});
                uint subdiv[3] = {0, i_vtx-1, i_vtx};

                // mappings from quadrature ref. el (bary) to intersection ref. elements (local coords)
                arma::mat master_map(master_dim, qdim + 1);
                arma::mat slave_map(slave_dim, qdim+1);
                for(uint i_col=0; i_col < qdim+1; i_col++) {
                    uint ip = subdiv[i_col];
                    master_map.col(i_col) = (*il)[ip].comp_coords();
                    slave_map.col(i_col) = (*il)[ip].bulk_coords();
                }

                double m_jac =
                    master_map.at(0,0)*(master_map.at(1,1) - master_map.at(1,2)) +
                    master_map.at(0,1)*(master_map.at(1,2) - master_map.at(1,0)) +
                    master_map.at(0,2)*(master_map.at(1,0) - master_map.at(1,1));

                double master_map_jac = arma::det( master_map.cols(1,2) - master_map.col(0) * arma::ones(1,2) );
                //DebugOut().fmt("jac1: {} jac2: {}", m_jac, master_map_jac);

                QGauss<qdim> q(2);
                arma::vec values(n_dofs);

                double q_area = 0;
                for(uint ip=0; ip < q.size(); ip++) {
                    arma::vec bary_qp =  RefElement<qdim>::local_to_bary(q.point(ip));

                    arma::vec qp_master = RefElement<master_dim>::local_to_bary(master_map * bary_qp);
                    arma::vec qp_slave = RefElement<slave_dim>::local_to_bary(slave_map * bary_qp);
                    double JxW = q.weight(ip) * isec_measure * master_jac;
                    uint shift = 0;
                    values.rows(shift, shift + master_dim) = P1_for_P1e[master_dim] * qp_master;
                    shift=master_dim + 1;
                    values.rows(shift, shift + slave_dim ) = -P1_for_P1e[slave_dim] * qp_slave;

                    double scale = -master_sigma * JxW;
                    arma::mat x = scale * ( values * values.t() );

                    //DebugOut() << print_var(scale);
                    //DebugOut() << print_var(values);
                    //DebugOut().fmt("2d one: {} 3d one: {}", arma::sum(values.rows(0,2)), arma::sum(values.rows(3,6)) );
                    //DebugOut() << x;
                    loc_system_.get_matrix() += x;

                    isec_area += q.weight(ip) * master_map_jac;
                    q_area += q.weight(ip);
                }


                //DebugOut().fmt("isub: {} isec m: {} a: {} qa: {} jac: {}\n",
                //        i_vtx, isec_measure, isec_area, q_area, master_map_jac);

            }

            loc_system_.eliminate_solution();
            if (fix_velocity_flag) {
                //this->fix_velocity_local(row_ele, col_ele);
            } else {
                data_->lin_sys->set_local_system(loc_system_);
            }

        } else {
            ASSERT(false).error("Impossible case.");
        }


    }
}


PL_CouplingAssembler::PL_CouplingAssembler(AssemblyDataPtr data)
: MortarAssemblyBase(data),
  slave_ac_(data->mh_dh),
  master_ac_(data->mh_dh)
{
    /// TODO: Use RefElement numberings to make it independent.
    P1_for_P1e[0] = arma::mat({{1}});
    P1_for_P1e[1] = arma::eye(2,2);

    //{{1, 1, -1}, {1, -1, 1}, {-1, 1, 1}}
    P1_for_P1e[2] = arma::ones(3,3);
    P1_for_P1e[2].at(0,2)=-1;
    P1_for_P1e[2].at(1,1)=-1;
    P1_for_P1e[2].at(2,0)=-1;

    //{{1, 1, 1, -1}, {1, 1, -1, 1}, {1, -1, 1, 1}, {-1, 1, 1, 1}}
    P1_for_P1e[3] = arma::ones(4,4);
    P1_for_P1e[3].at(0,3)=-2;
    P1_for_P1e[3].at(1,2)=-2;
    P1_for_P1e[3].at(2,1)=-2;
    P1_for_P1e[3].at(3,0)=-2;

}


void PL_CouplingAssembler::set_side(LocalElementAccessorBase<3> ele_ac,
        uint i_side, uint shift, arma::uvec &dofs)
{
    uint loc_dof = shift+i_side;
    dofs[loc_dof] =  ele_ac.edge_row(i_side);
    Boundary * bcd = ele_ac.full_iter()->side(i_side)->cond();

    if (bcd) {
        ElementAccessor<3> b_ele = bcd->element_accessor();
        auto type = (DarcyMH::EqData::BC_Type)data_->bc_type.value(b_ele.centre(), b_ele);
        //DebugOut().fmt("bcd id: {} sidx: {} type: {}\n", ele->id(), i_side, type);
        if (type == DarcyMH::EqData::dirichlet) {
            double bc_pressure = data_->bc_pressure.value(b_ele.centre(), b_ele);
            loc_system_.set_solution(loc_dof, bc_pressure, -1.0);
        }
    }
}

/*
template<uint qdim, uint master_dim, uint slave_dim>
void PL_CouplingAssembler::isec_assembly(
        double master_sigma,
        const IntersectionLocal<master_dim, slave_dim> &il,
        std::array<uint, qdim+1 > subdiv)
{
    // mappings from quadrature ref. el (bary) to intersection ref. elements (local coords)
    arma::mat master_map(master_dim, qdim + 1);
    arma::mat slave_map(slave_dim, qdim+1);
    for(uint i_col=0; i_col < qdim+1; i_col++) {
        uint ip = subdiv[i_col];
        master_map.col(i_col) = il[ip].comp_coords();
        slave_map.col(i_col) = il[ip].bulk_coords();
    }


    //QGauss<qdim> q(2);
    //arma::vec values(n_dofs);

    for(uint ip=0; ip < q.size(); ip++) {
        arma::vec bary_qp =  RefElement<qdim>::local_to_bary(q.point(ip));

        arma::vec qp_master = RefElement<master_dim>::local_to_bary(master_map * bary_qp);
        arma::vec qp_slave = RefElement<slave_dim>::local_to_bary(slave_map * bary_qp);
        double JxW = q.weight(ip) * il.compute_measure() * master_jac;
        uint shift = 0;
        values.rows(shift, shift + master_dim) = P1_for_P1e[master_dim] * qp_master;
        shift=master_dim + 1;
        values.rows(shift, shift + slave_dim ) = -P1_for_P1e[slave_dim] * qp_slave;

        double scale = -master_sigma * JxW;
        arma::mat x = scale * ( values * values.t() );

        //DebugOut() << print_var(scale);
        //DebugOut() << print_var(values);
        //DebugOut() << x;
        loc_system_.get_matrix() += x;

    }

    loc_system_.eliminate_solution();
    if (fix_velocity_flag) {
        //this->fix_velocity_local(row_ele, col_ele);
    } else {
        data_->lin_sys->set_local_system(loc_system_);
    }


}
*/


/**
 * P1 connection of different dimensions
 *
 * - 20.11. 2014 - very poor convergence, big error in pressure even at internal part of the fracture
 */

void PL_CouplingAssembler::assembly(LocalElementAccessorBase<3> master_ac) {

    if (master_ac.dim() > 2) return; // supported only for 1D and 2D master elements
    auto &isec_list = mixed_mesh_.element_intersections_[master_ac.ele_global_idx()];
    if (isec_list.size() == 0) return; // skip empty masters

    //slave_ac_.setup(master_ac);

    ElementFullIter master_ele = master_ac.full_iter();
    master_ac_.reinit(master_ac.ele_global_idx());
    arma::vec3 ele_centre = master_ele->centre();
    double m_sigma = data_->sigma.value( ele_centre, master_ele->element_accessor());
    double m_conductivity = data_->conductivity.value( ele_centre, master_ele->element_accessor());
    double m_crossection = data_->cross_section.value( ele_centre, master_ele->element_accessor() );

    double master_sigma = 2*m_sigma*m_conductivity *
                    2/ m_crossection;
                    /**
                     * ?? How to deal with anisotropy ??
                     * 3d-2d : compute nv of 2d triangle
                     * 2d-2d : interpret as 2d-1d-2d, should be symmetric master-slave
                     * 2d-1d : nv is tangent to 2d and normal to 1d
                    arma::dot(data_->anisotropy.value( ele_centre, ele->element_accessor())*nv, nv)
                    */

    double master_measure = master_ele->measure();
    DebugOut().fmt("master ele: {}  ", master_ac.ele_global_idx());

    for(auto &isec_pair : isec_list) {
        IntersectionLocalBase *isec = isec_pair.second;
        ASSERT_EQ_DBG(master_ac.ele_global_idx(), isec->component_ele_idx());

        //ElementFullIter slave_ele = data_->mesh->element(isec->bulk_ele_idx());
        slave_ac_.reinit(isec->bulk_ele_idx());

        if (typeid(*isec) == typeid(IntersectionLocal<1,2>)) {
            const uint m_dim = 1;
            const uint s_dim = 2;

            auto il = static_cast<const IntersectionLocal<1,2> *>(isec);
            ASSERT_EQ_DBG( il->size(), 2);
            ASSERT_EQ_DBG( master_ac_.dim(), 1);
            DebugOut() << "il: " <<  (*il)[0].comp_coords() << " " << (*il)[1].comp_coords();
            DebugOut() << "il: " <<  (*il)[0].bulk_coords() << " " << (*il)[1].bulk_coords();

            arma::mat::fixed<3, s_dim+1> ref_map_slave = elm_map2.element_map( *slave_ac_.full_iter() );
            arma::mat::fixed<3, m_dim+1> ref_map_master = elm_map1.element_map( *master_ac_.full_iter() );
            for(uint i_side=0; i_side< master_ac_.n_sides(); i_side ++) {

                auto side_center_ref = RefElement<m_dim>::local_to_bary(RefElement<m_dim>::centers_of_subelements(0)[i_side]);
                arma::vec3 side_center = ref_map_master * side_center_ref;
                arma::Col<double>::fixed<3>  bary_slave = elm_map2.project_real_to_unit( side_center, ref_map_slave);
                DebugOut() << "msc: " << side_center_ref << "b_slave:" << bary_slave;
                if (bary_slave.min() < -1e-6) continue;

                DebugOut().fmt("is: {} ",  i_side);

                // side center in slave element

                uint n_dofs = 1 + slave_ac_.n_sides();
                arma::uvec dofs(n_dofs);
                loc_system_.reset(n_dofs, n_dofs);
                set_side(master_ac_, i_side, 0, dofs);
                for(uint i_s_side=0; i_s_side < slave_ac_.n_sides(); i_s_side++)
                    set_side(slave_ac_, i_s_side, 1, dofs);
                loc_system_.set_dofs(dofs, dofs);

                arma::vec values = arma::ones(n_dofs, 1);
                values.rows(1,3) = -P1_for_P1e[2] * bary_slave;
                loc_system_.get_matrix() = -(master_sigma * master_measure / master_ac_.dim())
                        * ( values * values.t() );

                loc_system_.eliminate_solution();

                if (fix_velocity_flag) {
                    //this->fix_velocity_local(row_ele, col_ele);
                } else {
                    data_->lin_sys->set_local_system(loc_system_);
                }

            }
        } else
        if (typeid(*isec) == typeid(IntersectionLocal<2,2>)) {
            auto il = static_cast<const IntersectionLocal<2,2> *>(isec);
            ASSERT_EQ_DBG( il->size(), 2);
            ASSERT(false).error("22 not implemented");
            DebugOut() << "il22";

            //this->isec_assembly<1,2,2>(master_sigma, *il, {0,1});

        } else if (typeid(*isec) == typeid(IntersectionLocal<2,3>)) {
            auto il = static_cast<const IntersectionLocal<2,3> *>(isec);
            if (il->size() <= 2) continue; // skip degenerated intersections
            ASSERT(false).error("23 not implemented");
        } else {
            ASSERT(false).error("Impossible case.");
        }

    }
}







P01_CouplingAssembler::P01_CouplingAssembler(AssemblyDataPtr data)
: MortarAssemblyBase(data),
  tensor_average_(16),
  col_average_(4),
  quadrature_(*(data->mesh)),
  slave_ac_(data->mh_dh)
{
    isec_data_list.reserve(30);

    /// TODO: Use RefElement numberings to make it independent.
    P1_for_P1e[0] = arma::mat({{1}});
    P1_for_P1e[1] = arma::eye(2,2);

    //{{1, 1, -1}, {1, -1, 1}, {-1, 1, 1}}
    P1_for_P1e[2] = arma::ones(3,3);
    P1_for_P1e[2].at(0,2)=-1;
    P1_for_P1e[2].at(1,1)=-1;
    P1_for_P1e[2].at(2,0)=-1;

    //{{1, 1, 1, -2}, {1, 1, -2, 1}, {1, -2, 1, 1}, {-2, 1, 1, 1}}
    P1_for_P1e[3] = arma::ones(4,4);
    P1_for_P1e[3].at(0,3)=-2;
    P1_for_P1e[3].at(1,2)=-2;
    P1_for_P1e[3].at(2,1)=-2;
    P1_for_P1e[3].at(3,0)=-2;
}


void P01_CouplingAssembler::add_intersection(LocalElementAccessorBase<3> ele_ac) {
    isec_data_list.push_back(IsecData());
    IsecData &i_data = isec_data_list.back();

    i_data.dim= ele_ac.dim();
    //i_data.delta = 0;
    i_data.dofs.resize(ele_ac.n_sides());
    //i_data.vel_dofs.resize(ele_ac.n_sides());
    i_data.dirichlet_dofs.resize(ele_ac.n_sides());
    i_data.dirichlet_sol.resize(ele_ac.n_sides());
    i_data.n_dirichlet=0;
    i_data.z_values_.resize(ele_ac.n_sides());
    i_data.values_.zeros(slave_ac_.n_sides(),1);

    //i_data.ele_z_coord_=ele_ac.centre()[2];

    for(unsigned int i_side=0; i_side < ele_ac.n_sides(); i_side++ ) {
        i_data.dofs[i_side]=ele_ac.edge_row(i_side);
        i_data.z_values_[i_side] = arma::dot(data_->gravity_vec_, ele_ac.side(i_side)->centre());
        //i_data.vel_dofs[i_side] = ele_ac.side_row(i_side);
        //i_data.z_sides[i_side]=ele_ac.side(i_side)->centre()[2];
        //DebugOut().fmt("edge: {} {}", i_side, ele_ac.edge_row(i_side));
        Boundary * bcd = ele_ac.full_iter()->side(i_side)->cond();
        if (bcd) {
            ElementAccessor<3> b_ele = bcd->element_accessor();
            auto type = (DarcyMH::EqData::BC_Type)data_->bc_type.value(b_ele.centre(), b_ele);
            //DebugOut().fmt("bcd id: {} sidx: {} type: {}\n", ele->id(), i_side, type);
            if (type == DarcyMH::EqData::dirichlet) {
                //DebugOut().fmt("Dirichlet: {}\n", ele->index());
                double bc_pressure = data_->bc_pressure.value(b_ele.centre(), b_ele);
                i_data.dirichlet_dofs[i_data.n_dirichlet] = i_side;
                i_data.dirichlet_sol[i_data.n_dirichlet] = bc_pressure;
                i_data.n_dirichlet++;
            }
        }
    }
    //DebugOut().fmt("s: {} z: {}", i_data.dofs.n_rows, i_data.z_values_.n_rows);

}


template<uint qdim>
double jacobian(arma::mat master_map) {
    return arma::norm(master_map.col(1)-master_map.col(0), 2);
}

template<>
double jacobian<2>(arma::mat master_map) {
    return  fabs(master_map.at(0,0)*(master_map.at(1,1) - master_map.at(1,2)) +
            master_map.at(0,1)*(master_map.at(1,2) - master_map.at(1,0)) +
            master_map.at(0,2)*(master_map.at(1,0) - master_map.at(1,1)));

}



template<uint qdim, uint master_dim, uint slave_dim>
void P01_CouplingAssembler::integrate(const IntersectionLocal<master_dim, slave_dim> &il, std::array<uint, qdim+1 > subdiv)
{
    if (il.size() < 2) return   ;
    IsecData & isec_data = isec_data_list.back();

    // local coordinates of simplex integration domain vertices (columns)
    arma::mat slave_map(slave_dim, qdim+1);
    arma::mat master_map(master_dim, qdim+1);
    for(uint i_col=0; i_col < qdim+1; i_col++) {
        uint ip = subdiv[i_col];
        ASSERT_LT(ip, il.size());
        slave_map.col(i_col) = il[ip].bulk_coords();
        master_map.col(i_col) = il[ip].comp_coords();
    }

    QGauss<qdim> q(1);

    isec_data.delta = il.compute_measure();    // TODO: remove
    for(uint iq=0; iq < q.size(); iq++) {
        arma::vec bary_qp =  RefElement<qdim>::local_to_bary(q.point(iq));
        // barycentric coords of QP on slave
        arma::vec qp_slave = RefElement<slave_dim>::local_to_bary(slave_map * bary_qp);

        // weights for quadrature on subset of slave reference element

        double JxW = q.weight(iq) * jacobian<qdim>(master_map)
                     * double(qdim);    // We integrate over triangle, normalize to one.
        isec_sum += JxW;

        // part of integration over slave
        isec_data.values_ -= P1_for_P1e[slave_dim] * qp_slave * JxW;

        //DebugOut() << print_var(JxW) << print_var(isec_data.delta);
        //DebugOut() << print_var(qp_slave);
        //DebugOut() << x;

    }

}

//void P01_CouplingAssembler::add_intersection(LocalElementAccessorBase<3> ele_ac, double delta) {
//}




/**
 * Works well but there is large error next to the boundary.
 */
void P01_CouplingAssembler::assembly(LocalElementAccessorBase<3> master_ac)
{


    // on the intersection element we consider
    // * intersection dofs for master and slave
    //   those are dofs of the space into which we interpolate
    //   base functions from individual master and slave elements
    //   For the master dofs both are usualy eqivalent.
    // * original dofs - for master same as intersection dofs, for slave
    //   all dofs of slave elements

    // form list of intersection dofs, in this case pressures in barycenters
    // but we do not use those form MH system in order to allow second schur somplement. We rather map these
    // dofs to pressure traces, i.e. we use averages of traces as barycentric values.

    /**
     *  TODO:
     *  - pass through the master and all slaves and collect global dofs , bcd, solution.
     *    I.e. call Nx pressure_diff not NxNx.
     *
     *  - Is it safe to have duplicate rows in local_system?
     *  - Is it better to have more smaller local system then single big one?
     *  - use one big or more smaller local systems to set.
     */


    if (master_ac.dim() > 2) return; // supported only for 1D and 2D master elements
    auto &isec_list = mixed_mesh_.element_intersections_[master_ac.ele_global_idx()];
    if (isec_list.size() == 0) return; // skip empty masters

    //slave_ac_.setup(master_ac);

    ElementFullIter ele = master_ac.full_iter();
    arma::vec3 ele_centre = ele->centre();
    double m_sigma = data_->sigma.value( ele_centre, ele->element_accessor());
    double m_conductivity = data_->conductivity.value( ele_centre, ele->element_accessor());
    double m_crossection = data_->cross_section.value( ele_centre, ele->element_accessor() );

    double master_sigma = 2*m_sigma*m_conductivity *
                    2/ m_crossection;
                    /**
                     * ?? How to deal with anisotropy ??
                     * 3d-2d : compute nv of 2d triangle
                     * 2d-2d : interpret as 2d-1d-2d, should be symmetric master-slave
                     * 2d-1d : nv is tangent to 2d and normal to 1d
                    arma::dot(data_->anisotropy.value( ele_centre, ele->element_accessor())*nv, nv)
                    */





    // 2d-3d, 1d-2d

    isec_data_list.clear();
    //double cs_sqr_avg = 0.0;
    isec_sum = 0.0;
    for(auto &isec_pair : isec_list) {
        //double cs = data_->cross_section.value(slave_ac_.full_iter()->centre(), slave_ac_.full_iter()->element_accessor());

        IntersectionLocalBase *isec = isec_pair.second;
        ASSERT_EQ_DBG(master_ac.ele_global_idx(), isec->component_ele_idx());

        //ElementFullIter slave_ele = data_->mesh->element(isec->bulk_ele_idx());
        slave_ac_.reinit(isec->bulk_ele_idx());
        add_intersection(slave_ac_);
        if (typeid(*isec) == typeid(IntersectionLocal<1,2>)) {
            const IntersectionLocal<1,2> &il = *(static_cast<const IntersectionLocal<1,2> *>(isec));
            //ASSERT_EQ_DBG( il.size(), 2)(il.component_ele_idx())(il.bulk_ele_idx())(il);
            ASSERT_EQ(slave_ac_.dim(), 2);
            //DebugOut() << "il size: " << il.size();
            integrate<1,1,2>(il, {0,1});

        } else  if (typeid(*isec) == typeid(IntersectionLocal<2,2>)) {
            auto il = static_cast<const IntersectionLocal<2,2> *>(isec);
            //ASSERT_EQ_DBG( il->size(), 2);
            //integrate<1,2,2>(il, {0,1});
            DebugOut() << "il22, not implemented";

            //this->isec_assembly<1,2,2>(master_sigma, *il, {0,1});

        } else if (typeid(*isec) == typeid(IntersectionLocal<2,3>)) {
            const IntersectionLocal<2,3> &il = *(static_cast<const IntersectionLocal<2,3> *>(isec));
            ASSERT_EQ(slave_ac_.dim(), 3);

            for(uint i_vtx=2; i_vtx< il.size();  i_vtx++) {
                //this->isec_assembly<2,2,3>(master_sigma, *il, {0, i_vtx1, i_vtx1+1});
                integrate<2,2,3>(il, {0, i_vtx-1, i_vtx});
            }
        }
    } // loop over intersections
    element_isec=isec_sum;

    //if ( ! (slave_ac_.dim() == 2 && master_ac.dim() ==2 ) ) {
        //DebugOut() << print_var(isec_sum);
        if (fabs(isec_sum - 1.0) > 1E-5) {
            double diff = (isec_sum - 1.0);
            string out;
            for(auto & isec : isec_list) {
                slave_ac_.reinit(isec.second->bulk_ele_idx());
                out += fmt::format(" {}", slave_ac_.full_iter()->id());
            }
            WarningOut() << print_var(diff)
                    << print_var(ele->id())
                    << endl
                    << out;

        }
    //}

    // compute master values
    add_intersection(master_ac);
    arma::vec row_avg = arma::vec(master_ac.dim()+1);
    row_avg.fill(isec_sum / (master_ac.dim()+1));
    isec_data_list.back().values_ = row_avg;


    //DebugOut().fmt( "cs2: {} d0: {}", cs_sqr_avg, delta_0);
   // master_sigma = master_sigma * (cs_sqr_avg / isec_sum)
   //         / isec_sum;


    //add_to_linsys(master_sigma * ele->measure()*0.25);
    add_to_linsys(master_sigma * ele->measure()/isec_sum/isec_sum);


    // 2d-2d
    //DebugOut() << "2d-2d";
    /*
    if (i < isec_list.size()) {
        isec_data_list.clear();
        isec_sum = 0.0;
        for(; i < isec_list.size(); ++i) {
                quadrature_.reinit(isec_list[i].second);
                slave_ac_.reinit( quadrature_.slave_idx() );
                double isec_measure = quadrature_.measure();
                isec_sum += isec_measure;
                //DebugOut().fmt("Assembly22: {} {} {}", ele->id(), slave_ac_.full_iter()->id(), isec_measure);
                pressure_diff(slave_ac_, isec_measure);
        }
        pressure_diff(master_ac, -isec_sum);

        master_sigma = 100 * m_conductivity/ m_crossection / isec_sum;

        add_to_linsys(master_sigma);
    }
    */
}

void P01_CouplingAssembler::add_to_linsys(double scale)
{
    //DebugOut() << "New list" << print_var(isec_data_list.size());
    // rows
     for(IsecData &row_ele : isec_data_list) {

         //DebugOut() << row_ele.values_;
         //DebugOut() << row_ele.z_values_;
         //columns
         for(IsecData &col_ele : isec_data_list) {
             //DebugOut().fmt("s: {} z: {}", col_ele.dofs.n_rows, col_ele.z_values_.n_rows);
             arma::mat loc_mat =  -scale * row_ele.values_ * col_ele.values_.t();  // tensor product
             loc_system_.reset(row_ele.dofs, col_ele.dofs);

             for(uint i=0; i< row_ele.n_dirichlet; i++)
                 loc_system_.set_solution_row(row_ele.dirichlet_dofs[i], row_ele.dirichlet_sol[i], -1.0);
             for(uint i=0; i< col_ele.n_dirichlet; i++) loc_system_.set_solution_col(col_ele.dirichlet_dofs[i], col_ele.dirichlet_sol[i]);
             //ASSERT( arma::norm(product_,2) == 0.0 );
             loc_system_.set_matrix(loc_mat);
             // Must have z-coords for every side, can not use averaging vector

             arma::vec z_rhs = (-scale*arma::dot(col_ele.values_, col_ele.z_values_))*row_ele.values_;
             loc_system_.set_rhs( z_rhs );
             loc_system_.eliminate_solution();

             if (fix_velocity_flag) {
                 this->fix_velocity_local(row_ele, col_ele);
             } else {
                 data_->lin_sys->set_local_system(loc_system_);
             }
         }
     }
}


void P01_CouplingAssembler::fix_velocity_local(const IsecData &row_ele, const IsecData & col_ele)
{
    uint n_rows = row_ele.vel_dofs.n_rows;
    uint n_cols = col_ele.dofs.n_rows;
    arma::vec pressure(n_cols);
    arma::vec add_velocity(n_rows);
    double * solution = data_->lin_sys->get_solution_array();
    for(uint icol=0; icol < n_cols; icol++ ) pressure[icol] = solution[col_ele.dofs[icol]];
    add_velocity =  loc_system_.get_matrix() * pressure - loc_system_.get_rhs();
    DebugOut() << "fix_velocity\n" << pressure << add_velocity;
    for(uint irow=0; irow < n_rows; irow++ ) solution[row_ele.vel_dofs[irow]] += add_velocity[irow] ;
}
