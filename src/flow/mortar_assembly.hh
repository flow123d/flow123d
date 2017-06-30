/*
 * mortar_assembly.hh
 *
 *  Created on: Feb 22, 2017
 *      Author: jb
 */

#ifndef SRC_FLOW_MORTAR_ASSEMBLY_HH_
#define SRC_FLOW_MORTAR_ASSEMBLY_HH_

#include "mesh/mesh.h"
#include "quadrature/intersection_quadrature.hh"
#include "fem/mapping_p1.hh"
#include "flow/darcy_flow_mh.hh"
#include "la/local_system.hh"
#include <vector>


typedef std::shared_ptr<DarcyMH::EqData> AssemblyDataPtr;


class MortarAssemblyBase {
public:
    typedef vector<unsigned int> IsecList;

    MortarAssemblyBase(AssemblyDataPtr data)
    : data_(data),
      mixed_mesh_(data->mesh->mixed_intersections()),
      fix_velocity_flag(false),
      total_isec(0.0),
      element_isec(0.0)
    {

    }

    virtual ~MortarAssemblyBase() {
        //DebugOut() << "total isec: " << total_isec;
    };

    // Default assembly is empty to allow dummy implementation for dimensions without coupling.
    virtual void assembly(LocalElementAccessorBase<3> ele_ac) {};

    void fix_velocity(LocalElementAccessorBase<3> ele_ac) {
        fix_velocity_flag = true;
        this->assembly(ele_ac);
        total_isec+=element_isec;
        fix_velocity_flag = false;
    }

protected:
    AssemblyDataPtr data_;
    MixedMeshIntersections &mixed_mesh_;
    LocalSystem loc_system_;
    bool fix_velocity_flag;

    double total_isec;
    double element_isec;

};


struct IsecData {
    arma::uvec vel_dofs;
    arma::uvec dofs;
    unsigned int dim;
    double delta;
    double ele_z_coord_;
    arma::uvec dirichlet_dofs;
    arma::vec dirichlet_sol;
    unsigned int n_dirichlet;

    arma::vec values_;      // weights for side pressures
    arma::vec z_values_;    // Z values on sides
};


class P0_CouplingAssembler :public MortarAssemblyBase {
public:
    P0_CouplingAssembler(AssemblyDataPtr data);
    void assembly(LocalElementAccessorBase<3> ele_ac) override;
private:
    void fix_velocity_local(const IsecData & row_ele, const IsecData &col_ele);
    void pressure_diff(LocalElementAccessorBase<3> ele_ac, double delta);
    inline arma::mat & tensor_average(unsigned int row_dim, unsigned int col_dim) {
        return tensor_average_[4*row_dim + col_dim];
    }

    void add_to_linsys(double scale);

    vector<IsecData> isec_data_list;

    /// Row matrices to compute element pressure as average of boundary pressures
    std::vector< arma::mat > tensor_average_;
    std::vector< arma::vec > col_average_;
    IntersectionQuadratureP0 quadrature_;
    arma::mat product_;
    LocalElementAccessorBase<3> slave_ac_;


};



class P1_CouplingAssembler :public MortarAssemblyBase {
public:
    P1_CouplingAssembler(AssemblyDataPtr data);
    void assembly(LocalElementAccessorBase<3> ele_ac) override;
private:
    void add_sides(LocalElementAccessorBase<3> ele_ac, unsigned int shift, arma::uvec &dofs);

    template<uint qdim, uint master_dim, uint slave_dim>
    inline void isec_assembly(
            double master_sigma,
            const IntersectionLocal<master_dim, slave_dim> &il,
            std::array<uint, qdim+1 > subdiv);

    LocalElementAccessorBase<3> slave_ac_;
    LocalElementAccessorBase<3> master_ac_;


    /**
     * For every element dimension a transition matrix from the P1e base to
     * the P1 base in the space of linear polynomials.
     * P1e base have support points in side centers,
     * P1 basis have support points in vertices (corresponds to P1_disc FE)
     */
    arma::mat P1_for_P1e[4];
    double master_jac;


};


/**
 * Coupling based on edge support points interaction.
 *
 * sigma * ( p_1d -
 */
class PL_CouplingAssembler :public MortarAssemblyBase {
public:
    PL_CouplingAssembler(AssemblyDataPtr data);
    void assembly(LocalElementAccessorBase<3> ele_ac) override;
private:
    void set_side(LocalElementAccessorBase<3> ele_ac, uint i_side, uint shift, arma::uvec &dofs);

    template<uint qdim, uint master_dim, uint slave_dim>
    inline void isec_assembly(
            double master_sigma,
            const IntersectionLocal<master_dim, slave_dim> &il,
            std::array<uint, qdim+1 > subdiv);

    LocalElementAccessorBase<3> slave_ac_;
    LocalElementAccessorBase<3> master_ac_;


    /**
     * For every element dimension a transition matrix from the P1e base to
     * the P1 base in the space of linear polynomials.
     * P1e base have support points in side centers,
     * P1 basis have support points in vertices (corresponds to P1_disc FE)
     */
    arma::mat P1_for_P1e[4];
    arma::uvec master_dofs;
    uint n_master_dirichlet;
    arma::uvec master_dirichlet_dofs;
    arma::vec master_dirichlet_sol;

    double master_jac;
    MappingP1<1,3> elm_map1;
    MappingP1<2,3> elm_map2;
    MappingP1<3,3> elm_map3;

    //
    arma::Mat<double>::fixed<3,2> side_slave_point_12;
    arma::Mat<double>::fixed<4,3> side_slave_point_23;
};



class P01_CouplingAssembler :public MortarAssemblyBase {
public:
    P01_CouplingAssembler(AssemblyDataPtr data);
    void assembly(LocalElementAccessorBase<3> ele_ac) override;
private:
    void fix_velocity_local(const IsecData & row_ele, const IsecData &col_ele);
    void add_intersection(LocalElementAccessorBase<3> ele_ac);
    inline arma::mat & tensor_average(unsigned int row_dim, unsigned int col_dim) {
        return tensor_average_[4*row_dim + col_dim];
    }

    template<uint qdim, uint master_dim, uint slave_dim>
    void integrate(const IntersectionLocal<master_dim, slave_dim> &il, std::array<uint, qdim+1 > subdiv);

    void add_to_linsys(double scale);

    vector<IsecData> isec_data_list;

    /// Row matrices to compute element pressure as average of boundary pressures
    std::vector< arma::mat > tensor_average_;
    std::vector< arma::vec > col_average_;
    IntersectionQuadratureP0 quadrature_;
    arma::mat product_;
    LocalElementAccessorBase<3> slave_ac_;

    /**
     * For every element dimension a transition matrix from the P1e base to
     * the P1 base in the space of linear polynomials.
     * P1e base have support points in side centers,
     * P1 basis have support points in vertices (corresponds to P1_disc FE)
     */
    arma::mat P1_for_P1e[4];
    double isec_sum;

};


#endif /* SRC_FLOW_MORTAR_ASSEMBLY_HH_ */
