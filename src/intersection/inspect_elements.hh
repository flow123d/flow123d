/*
 * inspectelements.h
 *
 *  Created on: 13.4.2014
 *      Author: viktor, pe, jb
 */
#ifndef INSPECT_ELEMENTS_H_
#define INSPECT_ELEMENTS_H_

#include "inspect_elements_algorithm.hh"

class Mesh; // forward declare

namespace computeintersection {

class InspectElements;
class InspectElementsAlgorithm22;
class IntersectionLocalBase;
template<unsigned int N, unsigned int M> class IntersectionLocal;


enum IntersectionType
{
    none = 0,
    d12 = 0x0001,
    d13 = 0x0002,
    d23 = 0x0004,
    d22 = 0x0008,
    all = 0xFFFF
};

    
class InspectElements
{   
public:
    /// Stores 1D-3D intersections.
    std::vector<IntersectionLocal<1,3>> intersection_storage13_;
    /// Stores 2D-3D intersections.
    std::vector<IntersectionLocal<2,3>> intersection_storage23_;
    /// Stores 2D-2D intersections.
    std::vector<IntersectionLocal<2,2>> intersection_storage22_;
    
    /// Maps between elements and their intersections.
    /// i.e.: 
    /// intersection_map_[element index][i].first = other element index
    /// intersection_map_[element index][i].second = pointer to the intersection object
    std::vector<std::vector<ILpair>> intersection_map_;
    
    InspectElements(Mesh *mesh);
    ~InspectElements();
    
    /// Calls @p InspectElementsAlgorithm<dim>, computes intersections, 
    /// move them to storage, create the map and throw away the rest.
    void compute_intersections(IntersectionType d = IntersectionType::all);
    
    //temporary functions:
    
    /// Computes the size of the intersection in @p dim dimenstions.
    double measure_13();
    double measure_23();
    
    /// Generates a mesh file of the given name, including the intersection.
    void print_mesh_to_file_13(std::string name);
    void print_mesh_to_file_23(std::string name);
    
private:
    /// Mesh pointer.
    Mesh * mesh;
    
    InspectElementsAlgorithm<1> algorithm13_;
    InspectElementsAlgorithm<2> algorithm23_;
    InspectElementsAlgorithm22 algorithm22_;
    
    /// Auxiliary function that calls InspectElementsAlgorithm<dim>.
    template<unsigned int dim> void compute_intersections(InspectElementsAlgorithm<dim> &iea,
                                                          std::vector<IntersectionLocal<dim,3>> &storage);
    void compute_intersections_22(std::vector<IntersectionLocal<2,2>> &storage);
};

    
} // END namespace


#endif // INSPECT_ELEMENTS_H_