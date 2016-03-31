/*
 * inspectelements.h
 *
 *  Created on: 13.4.2014
 *      Author: viktor
 */
#ifndef INSPECT_ELEMENTS_H_
#define INSPECT_ELEMENTS_H_

#include "system/system.hh"

#include "mesh/bounding_box.hh"
#include "mesh/mesh_types.hh"

#include "simplex.h"

#include <queue>

class Mesh; // forward declare

namespace computeintersection {

template<unsigned int N, unsigned int M> class IntersectionPointAux;
template<unsigned int N, unsigned int M> class IntersectionAux;

class InspectElements;
class IntersectionLocalBase;
template<unsigned int N, unsigned int M> class IntersectionLocal;
template<unsigned int N, unsigned int M> class IntersectionPoint;



/** @brief Class implements algorithm for dim-dimensional intersections with 3D elements.
 * 
 * @p dim-D elements that are continuously neighbouring are called component elements.
 * 3D elements are called bulk elements.
 * 
 * Implements the initialization routine, that finds the first candidate for intersection.
 * It uses bounding boxes to fastly resolve intersection candidates. 
 * The lookup is done using BIH tree algorithm.
 * 
 * Implements prolongation algorithm that recursively searches neighbouring elements for next intersection candidates.
 * The candidates -- neighbouring component elements and bulk elements -- are pushed into separate queues.
 * The bulk elements queue is emptied at first, then the component elements queue is popped.
 * 
 * The recuring prolongation algorithm is as follows:
 * Function @p prolongation_decide fills the queues.
 * A candidate pair of elements is popped out of a queue.
 * Function @p prolongate computes intersection for a candidate pair and calls @p prolongation_decide again.
 * This is done in an infinite cycle, until both queues are empty.
 */
template<unsigned int dim>
class InspectElementsAlgorithm{
public:
    /** @brief Auxiliary structure for prolongation process.
     * 
     * Contains index of component element, index of bulk element 
     * and its index in vector of bulk elements intersecting with the component element.
     */
    struct Prolongation{
        unsigned int component_elm_idx;
        unsigned int elm_3D_idx;
        unsigned int dictionary_idx;
    };
    
    
    InspectElementsAlgorithm(Mesh *_mesh);
    ~InspectElementsAlgorithm();
    
    /// Runs the core algorithm for computing dimD-3D intersection.
    void compute_intersections();

private:
    const unsigned int undefined_elm_idx_ = -1;
    
    /// Counter for intersection among elements.
    unsigned int n_intersections_;
    
    /// Prolongation queue in the component mesh.
    std::queue<Prolongation> component_queue_;
    /// Prolongation queue in the bulk mesh.
    std::queue<Prolongation> bulk_queue_;
    
    // Array of flags, which elements are computed
    std::vector<bool> closed_elements;
    std::vector<unsigned int> last_slave_for_3D_elements;
    
    /// Mesh pointer.
    Mesh *mesh;
    /// Elements bounding boxes.
    std::vector<BoundingBox> elements_bb;
    /// Bounding box of all 3D elements.
    BoundingBox mesh_3D_bb;
    
    /// Object representing a single component element.
    Simplex<dim> component_simplex;
    /// Object representing a single bulk element.
    Simplex<3> tetrahedron;
    
    /// Resulting vector of intersections.
    std::vector<std::vector<IntersectionAux<dim,3>>> intersection_list_;
    
    /// Auxiliary vector that is filled in tracing algorithm of polygons (in 2D)
    /// and then used in prolongation decision routines.
    std::vector<unsigned int> prolongation_table_;
    
    /// Initialization.
    /// Sets vector sizes and computes bulk bounding box.
    void init();
    
    /// Auxiliary function that translates @p ElementFullIter to @p Simplex<simplex_dim>.
    template<unsigned int simplex_dim>
    void update_simplex(const ElementFullIter &element, Simplex<simplex_dim> & simplex);
    
    /// A hard way to find whether the intersection of two elements has already been computed, or not.
    bool intersection_exists(unsigned int component_ele_idx, unsigned int bulk_ele_idx);

    /// Auxiliary function for calling tracing algorithm. Is empty if @p dim =0.
    void trace(IntersectionAux<dim,3> &intersection);
    
    /// Computes the first intersection, from which we then prolongate.
    bool compute_initial_CI(unsigned int component_ele_idx, unsigned int bulk_ele_idx);
    
    /// Finds neighbouring elements that are new candidates for intersection and pushes
    /// them into component queue or bulk queue.
    void prolongation_decide(const ElementFullIter &elm, const ElementFullIter &ele_3D);
    
    /// Computes the intersection for a candidate in a queue and calls @p prolongation_decide again.
    void prolongate(const Prolongation &pr);
    
    friend class InspectElements;
};




class InspectElements
{   
public:
    /// First = element index, Second = pointer to intersection object.
    typedef std::pair<unsigned int, IntersectionLocalBase*> ILpair;
    
    /// Stores 1D-3D intersections.
    std::vector<IntersectionLocal<1,3>> intersection_storage13_;
    /// Stores 2D-3D intersections.
    std::vector<IntersectionLocal<2,3>> intersection_storage23_;
    
    /// Maps between elements and their intersections.
    /// i.e.: 
    /// intersection_map_[element index][i].first = other element index
    /// intersection_map_[element index][i].second = pointer to the intersection object
    std::vector<std::vector<ILpair>> intersection_map_;
    
    InspectElements(Mesh *mesh);
    ~InspectElements();
    
    /// Calls @p InspectElementsAlgorithm<dim>, computes intersections, 
    /// move them to storage, create the map and throw away the rest.
    void compute_intersections();
    
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
    
    template<unsigned int dim> void compute_intersections(std::vector<IntersectionLocal<dim,3>> &storage);
};

    
} // END namespace


#endif // INSPECT_ELEMENTS_H_