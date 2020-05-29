/*
 * assembly_theory.hh
 *
 *  Created on: Apr 30, 2016
 *      Author: jb
 */

#ifndef SRC_FLOW_ASSEMBLY_THEORY_HH_
#define SRC_FLOW_ASSEMBLY_THEORY_HH_

/**
 * @file
 * Conceptual scheme of Assembly classes.
 *
 * TODO:
 *  - passing data between classes
 *  - equations derive code and data
 *  - assembly residual and assembly jacobian
 */


class AssemblyGroupInterface {

};

template <class ImplAL>
class AssemblyGroupBase : public AssemblyGroupInterface {
    void group_generate_quadratute_points();
    void group_precompute_fields();
    void group_assembly_element_integrals();
    void group_assembly_edge_integrals();
    void group_assembly_boundary_integrals();

    // interface for descendants


};




/**
 * Class that groups the mesh entities (mainly elements) into
 * local groups that should fit into the cache with all necessary data.
 * Assembly of one group may involve few virtual calls, and may involve more cycles, since
 * these are local in the cache.
 * Group assembly scheme:
 * 1. generate quadrature points, make subgroups by region (partly VECTORISE)
 * 2. precompute field values in quadrature points (may involve various dependencies)
 *    single virtual call per involved field per subgroup (VECTORISE)
 * 3. assembly element integrals per dim (VECTORISE)
 * 4. assembly edge integrals per dim (VECTORISE)
 * 5. assembly boundary integrals
 */
class Assembler {
    typedef std::vector< std::shared_ptr<AssemblyGroupInterface> > MultidimGroupAssembly;
    MultidimGroupAssembly group_assemblers_;
    Assembler(MultidimGroupAssembly assembly_vec_)
    : group_assemblers_(assembly_vec_)
    {}


};

/********************** Equations examples ********************/


class FlowMHData : public FieldSet {

};

template <int dim>
class FlowMHAssembly : public AssemblyGroupBase< SomeEqAssembly<dim> > {
    FlowMHAssembly(std::shared_ptr<FlowMHData> data);
    std::shared_ptr<FlowMHData> data_;
};

class FlowMHEq : public EquationBase {
    FlowMHEq(Mesh &mesh) {

    }
};

/********** Derived equation **********/

class FlowLMHData : public FlowMHData {

};

template <int dim>
class FlowLMHAssembly : public AssemblyGroupBase< SomeEqAssembly<dim> > {
    FlowLMHAssembly(std::shared_ptr<FlowLMHData> data)
    :FlowMHAssembly(data)
    {}

    void assembly_vector_on_element();

    void assembly_matrix_on_element();

    std::shared_ptr<FlowLMHData> data_;
};

class FlowLMHEq : public EquationBase {
    FlowLMHEq(Mesh &mesh) {

    }
};


#endif /* SRC_FLOW_ASSEMBLY_THEORY_HH_ */
