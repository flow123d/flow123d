/*
 * mortar_assembly.hh
 *
 *  Created on: Feb 22, 2017
 *      Author: jb
 */

#ifndef SRC_FLOW_MORTAR_ASSEMBLY_HH_
#define SRC_FLOW_MORTAR_ASSEMBLY_HH_

#include "mesh/mesh.h"
#include "flow/darcy_flow_mh.hh"


typedef std::shared_ptr<DarcyMH::EqData> AssemblyDataPtr;


class MortarAssemblyBase {
public:
    typedef vector<unsigned int> IsecList;

    // Default assembly is empty to allow dummy implementation for dimensions without coupling.
    virtual void assembly(LocalElementAccessorBase<3> ele_ac) {};
    MortarAssemblyBase(AssemblyDataPtr data)
    : data_(data),
      master_list_(data->mesh->master_elements),
      intersections_(data->mesh->intersections)
    {

    }

    virtual ~MortarAssemblyBase() {};

protected:
    AssemblyDataPtr data_;
    const vector<IsecList> &master_list_;
    const vector<Intersection> &intersections_;
    LocalSystem loc_system_;

};


struct IsecData {
    arma::uvec dofs;
    unsigned int dim;
    double delta;
    arma::uvec dirichlet_dofs;
    arma::vec dirichlet_sol;
    unsigned int n_dirichlet;
};


class P0_CouplingAssembler :public MortarAssemblyBase {
public:
    P0_CouplingAssembler(AssemblyDataPtr data);
    void assembly(LocalElementAccessorBase<3> ele_ac);
    void pressure_diff(uint i_ele, const IsecList &isec_list, LocalElementAccessorBase<3> ele_ac, IsecData &i_data);
private:
    inline arma::mat & tensor_average(unsigned int row_dim, unsigned int col_dim) {
        return tensor_average_[4*row_dim + col_dim];
    }

    vector<IsecData> isec_data_list;

    /// Row matrices to compute element pressure as average of boundary pressures
    vector< vector< arma::mat > > tensor_average;
    /// measure of master element, should be sum of intersection measures
    double delta_0;
};



class P1_CouplingAssembler :public MortarAssemblyBase {
public:
    P1_CouplingAssembler(AssemblyDataPtr data)
: MortarAssemblyBase(data),
      rhs(5),
      dofs(5),
      dirichlet(5)
    {
        rhs.zeros();
    }

    void assembly(LocalElementAccessorBase<3> ele_ac);
    void add_sides(LocalElementAccessorBase<3> ele_ac, unsigned int shift, vector<int> &dofs, vector<double> &dirichlet);
private:

    arma::vec rhs;
    vector<int> dofs;
    vector<double> dirichlet;
};



#endif /* SRC_FLOW_MORTAR_ASSEMBLY_HH_ */
