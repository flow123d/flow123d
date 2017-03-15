/*
 * mortar_assembly.cc
 *
 *  Created on: Feb 22, 2017
 *      Author: jb
 */

#include "flow/mortar_assembly.hh"
#include "quadrature/intersection_quadrature.hh"
#include "la/linsys.hh"
//#include "mesh/intersection.hh"
#include "intersection/mixed_mesh_intersections.hh"
#include <armadillo>

P0_CouplingAssembler::P0_CouplingAssembler(AssemblyDataPtr data)
: MortarAssemblyBase(data),
  tensor_average_(16),
  quadrature_(*(data->mesh))
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
        }
}


void P0_CouplingAssembler::pressure_diff(LocalElementAccessorBase<3> ele_ac, double delta) {
    isec_data_list.push_back(IsecData());
    IsecData &i_data = isec_data_list.back();

    i_data.dim= ele_ac.dim();
    i_data.delta = delta;
    i_data.dofs.resize(ele_ac.n_sides());
    i_data.dirichlet_dofs.resize(ele_ac.n_sides());
    i_data.dirichlet_sol.resize(ele_ac.n_sides());
    i_data.n_dirichlet=0;

    for(unsigned int i_side=0; i_side < ele_ac.n_sides(); i_side++ ) {
        i_data.dofs[i_side]=ele_ac.edge_row(i_side);
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
 void P0_CouplingAssembler::assembly(LocalElementAccessorBase<3> ele_ac) {
    arma::mat product;
    unsigned int i,j;

    if (ele_ac.dim() > 2) return; // supported only for 1D and 2D master elements

    auto &isec_list = mixed_mesh_.element_intersections_[ele_ac.ele_global_idx()];
    if (isec_list.size() == 0) return; // skip empty masters


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
    ElementFullIter ele = ele_ac.full_iter();
    arma::vec3 ele_centre = ele->centre();
    double master_sigma = 2*data_->sigma.value( ele_centre, ele->element_accessor()) *
                    2*data_->conductivity.value( ele_centre, ele->element_accessor())
                    /**
                     * ?? How to deal with anisotropy ??
                     * 3d-2d : compute nv of 2d triangle
                     * 2d-2d : interpret as 2d-1d-2d, should be symmetric master-slave
                     * 2d-1d : nv is tangent to 2d and normal to 1d
                    arma::dot(data_->anisotropy.value( ele_centre, ele->element_accessor())*nv, nv)
                    */
                    / data_->cross_section.value( ele_centre, ele->element_accessor() );       // crossection of lower dim.

                    // data_->cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) * // cross-section of higher dim. (2d)
                    // data_->cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) /




    uint master_dim = ele_ac.dim();
    uint m_idx = ele_ac.full_iter()->id();
    isec_data_list.clear();

    pressure_diff(ele_ac, -1.0);
    double cs_sqr_avg = 0.0;
    double isec_sum = 0.0;
    for(i = 0; i < isec_list.size(); ++i) {
        quadrature_.reinit(isec_list[i].second);
        ele_ac.reinit( quadrature_.slave_idx() );
        double cs = data_->cross_section.value(ele_ac.full_iter()->centre(), ele_ac.full_iter()->element_accessor());
        double isec_measure = quadrature_.measure();
        //DebugOut() << print_var(cs) << print_var(isec_measure);
        cs_sqr_avg += cs*cs*isec_measure;
        isec_sum += isec_measure;
        //DebugOut().fmt("Assembly: {} {} {}", m_idx, ele_ac.full_iter()->id(), quadrature_.measure());
        pressure_diff(ele_ac, isec_measure);
    }
    if ( ! (ele_ac.dim() == 2 && master_dim ==2 ) ) {
        double ref_master_measure = 1.0 / master_dim;
        if ( abs(isec_sum - ref_master_measure) > 1E-5)
            WarningOut().fmt("Wrong intersection area: {} != {}", isec_sum, ref_master_measure );
    }
    isec_data_list[0].delta = -isec_sum;

    //DebugOut().fmt( "cs2: {} d0: {}", cs_sqr_avg, delta_0);
    master_sigma = master_sigma * (cs_sqr_avg / isec_sum)
            / isec_sum;


    // rows
    for(IsecData &row_ele : isec_data_list) {
        //columns
        for(IsecData &col_ele : isec_data_list) {


            double scale =  -master_sigma * row_ele.delta * col_ele.delta;
            product = scale * tensor_average(row_ele.dim, col_ele.dim);

            loc_system_.reset(row_ele.dofs, col_ele.dofs);

            for(i=0; i< row_ele.n_dirichlet; i++) loc_system_.set_solution_row(row_ele.dirichlet_dofs[i], row_ele.dirichlet_sol[i]);
            for(i=0; i< col_ele.n_dirichlet; i++) loc_system_.set_solution_col(col_ele.dirichlet_dofs[i], col_ele.dirichlet_sol[i]);
            loc_system_.set_matrix(product);
            data_->lin_sys->set_local_system(loc_system_);
        }
    }
 }



 void P1_CouplingAssembler::add_sides(LocalElementAccessorBase<3> ele_ac, unsigned int shift, vector<int> &dofs, vector<double> &dirichlet)
 {
        for(unsigned int i_side=0; i_side < ele_ac.n_sides(); i_side++ ) {
            dofs[shift+i_side] =  ele_ac.edge_row(i_side);
            Boundary * bcd = ele_ac.full_iter()->side(i_side)->cond();

            if (bcd) {
                ElementAccessor<3> b_ele = bcd->element_accessor();
                auto type = (DarcyMH::EqData::BC_Type)data_->bc_type.value(b_ele.centre(), b_ele);
                //DebugOut().fmt("bcd id: {} sidx: {} type: {}\n", ele->id(), i_side, type);
                if (type == DarcyMH::EqData::dirichlet) {
                    //DebugOut().fmt("Dirichlet: {}\n", ele->index());
                    dofs[shift + i_side] = -dofs[shift + i_side];
                    double bc_pressure = data_->bc_pressure.value(b_ele.centre(), b_ele);
                    dirichlet[shift + i_side] = bc_pressure;
                }
            }
        }
 }


/**
 * P1 connection of different dimensions
 *
 * - 20.11. 2014 - very poor convergence, big error in pressure even at internal part of the fracture
 */

void P1_CouplingAssembler::assembly(LocalElementAccessorBase<3> ele_ac) {
/*
    const IsecList &master_list = master_list_[ele_ac.ele_global_idx()];
    if (master_list.size() == 0) return; // skip empty masters
    double master_sigma=data_->sigma.value( ele_ac.full_iter()->centre(), ele_ac.element_accessor());

    // set mater sides
    add_sides(ele_ac, 3, dofs, dirichlet);

    for(uint i = 0; i < master_list.size(); ++i) {
        const IntersectionQuadrature &intersec = intersections_[ master_list[i] ];
        const Element * slave = intersec.slave_iter();
        ele_ac.reinit(intersec.slave_iter()->index());
        add_sides(ele_ac, 0, dofs, dirichlet);
*/

/*
 * Local coordinates on 1D
 *         t0
 * node 0: 0.0
 * node 1: 1.0
 *
 * base fce points
 * t0 = 0.0    on side 0 node 0
 * t0 = 1.0    on side 1 node 1
 *
 * Local coordinates on 2D
 *         t0  t1
 * node 0: 0.0 0.0
 * node 1: 1.0 0.0
 * node 2: 0.0 1.0
 *
 * base fce points
 * t0=0.5, t1=0.0        on side 0 nodes (0,1)
 * t0=0.5, t1=0.5        on side 1 nodes (1,2)
 * t0=0.0, t1=0.5        on side 2 nodes (2,0)
 */

/*

        arma::vec point_Y(1);
        point_Y.fill(1.0);
        arma::vec point_2D_Y(intersec.map_to_slave(point_Y)); // local coordinates of  Y on slave (1, t0, t1)
        arma::vec point_1D_Y(intersec.map_to_master(point_Y)); //  local coordinates of  Y on master (1, t0)

        arma::vec point_X(1);
        point_X.fill(0.0);
        arma::vec point_2D_X(intersec.map_to_slave(point_X)); // local coordinates of  X on slave (1, t0, t1)
        arma::vec point_1D_X(intersec.map_to_master(point_X)); // local coordinates of  X on master (1, t0)

        arma::mat base_2D(3, 3);
        // basis functions are numbered as sides
        // TODO:
        // Use RT finite element to evaluate following matrices.

        // Ravirat - Thomas base functions evaluated in points (0,0), (1,0), (0,1)
        // 2D RT_i(t0, t1) = a0 + a1*t0 + a2*t1
        //         a0     a1      a2
        base_2D << 1.0 << 0.0 << -2.0 << arma::endr // RT for side 0
                << 1.0 << -2.0 << 0.0 << arma::endr // RT for side 1
                << -1.0 << 2.0 << 2.0 << arma::endr;// RT for side 2


        arma::mat base_1D(2, 2);
        // Ravirat - Thomas base functions evaluated in points (0,0), (1,0), (0,1)
        // 1D RT_i(t0) =   a0 + a1 * t0
        //          a0     a1
        base_1D << 1.0 << -1.0 << arma::endr // RT for side 0,
                << 0.0 << 1.0 << arma::endr; // RT for side 1,



        // Consider both 2D and 1D value are defined for the test function
        // related to the each of 5 DOFs involved in the intersection.
        // One of these values is always zero.
        // Compute difference of the 2D and 1D value for every DOF.
        // Compute value of this difference in both endpoints X,Y of the intersection.

        arma::vec difference_in_Y(5);
        arma::vec difference_in_X(5);
        // slave sides 0,1,2
        difference_in_Y.subvec(0, 2) = -base_2D * point_2D_Y;
        difference_in_X.subvec(0, 2) = -base_2D * point_2D_X;
        // master sides 3,4
        difference_in_Y.subvec(3, 4) = base_1D * point_1D_Y;
        difference_in_X.subvec(3, 4) = base_1D * point_1D_X;

        // applying the Simpson's rule
        // to the product of two linear functions f, g we get
        // (b-a)/6 * ( 3*f(a)*g(a) + 3*f(b)*g(b) + 2*f(a)*g(b) + 2*f(b)*g(a) )
        arma::mat A(5, 5);
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                A(i, j) = -master_sigma * intersec.intersection_true_size() *
                        ( difference_in_Y[i] * difference_in_Y[j]
                          + difference_in_Y[i] * difference_in_X[j]/2
                          + difference_in_X[i] * difference_in_Y[j]/2
                          + difference_in_X[i] * difference_in_X[j]
                        ) * (1.0 / 3.0);

            }
        }
        auto dofs_cp=dofs;
        data_->lin_sys->set_values( dofs_cp, dofs_cp, A, rhs, dirichlet, dirichlet);

    }*/
}



