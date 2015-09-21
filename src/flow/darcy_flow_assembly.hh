/*
 * darcy_flow_assembly.hh
 *
 *  Created on: Sep 21, 2015
 *      Author: jb
 */

#ifndef SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_
#define SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_

#include <type_traits>

namespace darcy_flow {

    /**
     * Prototype of an abstract local assembly class.
     * TODO: Make it independent of equation.
     */
    class LocalAssemblyBase
    {
    public:
        // assembly just A block of local matrix
        virtual void assembly_local_matrix(arma::mat &local_matrix,
                                           ElementFullIter ele) = 0;

        // assembly compatible neighbourings
        virtual void assembly_local_vb(double *local_vb,
                                       ElementFullIter ele,
                                       Neighbour *ngh) = 0;

        // compute velocity value in the barycenter
        // TOTO: implement and use general interpolations between discrete spaces
        virtual arma::vec3 barycenter_velocity(ElementFullIter ele) = 0;

        virtual ~LocalAssemblyBase();
    };

    /**
     * Prototype of a global assembly class. This takes care of iterating through the mesh
     * and performing possible cashing to optimize the assembly. Should contain only code
     * independent of particular equation.
     *
     * LocalAssembly template provides implementations of LocalAssemblyBase dependent on dimension.
     */
    template<template<typename> class LocalAssembly>
    class Assembler {
        typedef LocalAssembly<dim>::AssemblyData Data;
        Assembler(Data data) {
            set_data(data);
        }
        void set_data(Data data) {
            make_local_assembly<3>(data);
        }
        LocalAssemblyBase &assembly(unsigned int dim) {
            ASSERT( dim>0 && dim <4, "Wrong dimension.\n");
            return local_assembly[dim];
        }

    private:
        template<int dim>
        void make_local_assembly(Data data) {
            static_assert(std::is_base_of<LocalAssemblyBase, LocalAssembly<dim> >::value);
            local_assembly[dim] = new LocalAssembly<dim>(data);
            make_local_assembly<dim-1>(data);
        }
        template<>
        void make_local_assembly<0>(Data data) {}
        LocalAssemblyBase * local_assembly[3];
    };

    struct AssemblyData
    {
        Mesh *mesh;
        EqData* data;
        MH_DofHandler *mh_dh;
    };






} // namespace darcy_flow


#endif /* SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_ */
