/*
 * mortar_assembly.cc
 *
 *  Created on: Feb 22, 2017
 *      Author: jb
 */

#include "flow/mortar_assembly.hh"
#include "quadrature/intersection_quadrature.hh"
#include "la/linsys.hh"
#include "mesh/accessors.hh"
#include "fem/dh_cell_accessor.hh"
#include "intersection/mixed_mesh_intersections.hh"
#include <armadillo>

P0_CouplingAssembler::P0_CouplingAssembler(AssemblyDataPtr data)
: MortarAssemblyBase(data),
  tensor_average_(16),
  col_average_(4),
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
            col_average_[row_dim] = row_avg;
        }
}


void P0_CouplingAssembler::pressure_diff(const DHCellAccessor& dh_cell, double delta) {
    isec_data_list.push_back(IsecData());
    IsecData &i_data = isec_data_list.back();

    ElementAccessor<3> ele = dh_cell.elm();
    const uint nsides = ele->n_sides();
    const uint ndofs = dh_cell.n_dofs();

    i_data.dim= ele.dim();
    i_data.delta = delta;
    i_data.dofs.resize(nsides);
    i_data.vel_dofs.resize(nsides);
    i_data.dirichlet_dofs.resize(nsides);
    i_data.dirichlet_sol.resize(nsides);
    i_data.n_dirichlet=0;
    i_data.ele_z_coord_=ele.centre()[2];

    for(unsigned int i_side=0; i_side < nsides; i_side++ ) {
        // TODO: replace with DHCell getter when available for FESystem component
        i_data.dofs[i_side] = dh_cell.get_loc_dof_indices()[(ndofs+1)/2+i_side];   //edge dof
        i_data.vel_dofs[i_side] = dh_cell.get_loc_dof_indices()[i_side];   // side dof
        //i_data.z_sides[i_side]=ele.side(i_side)->centre()[2];
        Side side = *dh_cell.elm().side(i_side);
        if (side.is_boundary()) {
            Boundary bcd = side.cond();
            ElementAccessor<3> b_ele = bcd.element_accessor();
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
void P0_CouplingAssembler::assembly(const DHCellAccessor& dh_cell_master)
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

    ElementAccessor<3> ele = dh_cell_master.elm();
    
    if (ele.dim() > 2) return; // supported only for 1D and 2D master elements
    auto &isec_list = mixed_mesh_.element_intersections_[ele.idx()];
    if (isec_list.size() == 0) return; // skip empty masters

    //slave_ac_.setup(master_ac);

    arma::vec3 ele_centre = ele.centre();
    double m_sigma = data_->sigma.value( ele_centre, ele);
    double m_conductivity = data_->conductivity.value( ele_centre, ele);
    double m_crossection = data_->cross_section.value( ele_centre, ele );

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
    unsigned int slave_dim = 0;
    uint i = 0;
    for(; i < isec_list.size(); ++i) {
        bool non_zero = quadrature_.reinit(isec_list[i].second);
        DHCellAccessor dh_cell_slave(this->data_->dh_.get(), quadrature_.slave_idx());
        ElementAccessor<3> ele_slave = dh_cell_slave.elm();
        slave_dim = ele_slave.dim();
        if (slave_dim == ele.dim()) break;
        if (! non_zero) continue; // skip quadratures close to zero

        double cs = data_->cross_section.value(ele_slave.centre(), ele_slave);
        double isec_measure = quadrature_.measure();
        //DebugOut() << print_var(cs) << print_var(isec_measure);
        cs_sqr_avg += cs*cs*isec_measure;
        isec_sum += isec_measure;
        //DebugOut().fmt("Assembly23: {} {} {} ", ele.idx(), ele_slave->id(), isec_measure);
        pressure_diff(dh_cell_slave, isec_measure);
    }
    if ( ! (slave_dim == 2 && ele.dim() ==2 ) ) {
        if( fabs(isec_sum - ele.measure()) > 1E-5) {
            string out;
            for(auto & isec : isec_list) {
                DHCellAccessor dh_cell_slave(this->data_->dh_.get(), isec.second->bulk_ele_idx());
                out += fmt::format(" {}", dh_cell_slave.elm().idx());
            }

            double diff = (isec_sum - ele.measure())/ele.measure();
            WarningOut() << print_var(diff)
                    << print_var(ele.idx())
                    << endl
                    << out;

        }
    }
    pressure_diff(dh_cell_master, -isec_sum);

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
                DHCellAccessor dh_cell_slave(this->data_->dh_.get(), quadrature_.slave_idx());
                double isec_measure = quadrature_.measure();
                isec_sum += isec_measure;
                //DebugOut().fmt("Assembly22: {} {} {}", ele.idx(), dh_cell_slave.elm().idx(), isec_measure);
                pressure_diff(dh_cell_slave, isec_measure);
        }
        pressure_diff(dh_cell_master, -isec_sum);

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

             loc_system_.reset(row_ele.dofs, col_ele.dofs);

             for(uint i=0; i< row_ele.n_dirichlet; i++)
                 loc_system_.set_solution_row(row_ele.dirichlet_dofs[i], row_ele.dirichlet_sol[i], -1.0);
             for(uint i=0; i< col_ele.n_dirichlet; i++) loc_system_.set_solution_col(col_ele.dirichlet_dofs[i], col_ele.dirichlet_sol[i]);
             //ASSERT( arma::norm(product_,2) == 0.0 );
             loc_system_.set_matrix(product_);
             // Must have z-coords for every side, can not use averaging vector
             loc_system_.set_rhs( -s * col_average_[row_ele.dim] * col_ele.ele_z_coord_ );
             loc_system_.eliminate_solution();

             if (fix_velocity_flag) {
                 this->fix_velocity_local(row_ele, col_ele);
             } else {
                 loc_system_.eliminate_solution();
                 data_->lin_sys->set_local_system(loc_system_, data_->dh_->get_local_to_global_map());
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
    
    for(uint icol=0; icol < n_cols; icol++ ) pressure[icol] = data_->full_solution[col_ele.dofs[icol]];
    add_velocity =  loc_system_.get_matrix() * pressure - loc_system_.get_rhs();
    //DebugOut() << "fix_velocity\n" << pressure << add_velocity;
    for(uint irow=0; irow < n_rows; irow++ ) data_->full_solution[row_ele.vel_dofs[irow]] += add_velocity[irow] ;
}

void P1_CouplingAssembler::add_sides(const DHCellAccessor& dh_cell, unsigned int shift, vector<int> &dofs, vector<double> &dirichlet)
{
    ElementAccessor<3> ele = dh_cell.elm();
    const uint ndofs = dh_cell.n_dofs();

    for(unsigned int i_side=0; i_side < ele->n_sides(); i_side++ ) {
        // TODO: replace with DHCell getter when available for FESystem component
        dofs[shift+i_side] =  dh_cell.get_loc_dof_indices()[(ndofs+1)/2+i_side];   //edge dof
        
        Side side = *dh_cell.elm().side(i_side);
        if (side.is_boundary()) {
            Boundary bcd = side.cond();
            ElementAccessor<3> b_ele = bcd.element_accessor();
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

void P1_CouplingAssembler::assembly(FMT_UNUSED const DHCellAccessor& dh_cell) {
/*
    const IsecList &master_list = master_list_[ele_ac.ele_global_idx()];
    if (master_list.size() == 0) return; // skip empty masters
    double master_sigma=data_->sigma.value( ele_ac.element_accessor()->centre(), ele_ac.element_accessor());

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



