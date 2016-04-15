/*
 * inspectelements.h
 *
 *  Created on: 13.4.2014
 *      Author: viktor, pe, jb
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
class InspectElementsAlgorithm22;
class IntersectionLocalBase;
template<unsigned int N, unsigned int M> class IntersectionLocal;
template<unsigned int N, unsigned int M> class IntersectionPoint;

/// First = element index, Second = pointer to intersection object.
typedef std::pair<unsigned int, IntersectionLocalBase*> ILpair;

enum IntersectionType
{
    none = 0,
    d12 = 0x0001,
    d13 = 0x0002,
    d23 = 0x0004,
    d22 = 0x0008,
    all = 0xFFFF
};


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
    
    /// Initialization.
    /// Sets vector sizes and computes bulk bounding box.
    void init();
    
    /// Auxiliary function that translates @p ElementFullIter to @p Simplex<simplex_dim>.
    template<unsigned int simplex_dim>
    void update_simplex(const ElementFullIter &element, Simplex<simplex_dim> & simplex);
    
    /// A hard way to find whether the intersection of two elements has already been computed, or not.
    bool intersection_exists(unsigned int component_ele_idx, unsigned int bulk_ele_idx);
    
    /// Computes the first intersection, from which we then prolongate.
    bool compute_initial_CI(unsigned int component_ele_idx, unsigned int bulk_ele_idx, 
                            std::vector<unsigned int> &prolongation_table);
    
    /// Finds neighbouring elements that are new candidates for intersection and pushes
    /// them into component queue or bulk queue.
    void prolongation_decide(const ElementFullIter &comp_ele, const ElementFullIter &bulk_ele, 
                             const IntersectionAux<dim,3> &is,
                             const std::vector<unsigned int> &prolongation_table);
    
    /// Computes the intersection for a candidate in a queue and calls @p prolongation_decide again.
    void prolongate(const Prolongation &pr);
    
    std::vector<unsigned int> get_bulk_element_edges(const ElementFullIter& bulk_ele, 
                                                     const IntersectionPointAux<dim,3> &IP,
                                                     const bool &include_current_bulk_ele = false
                                                    );
    unsigned int create_prolongations_over_bulk_element_edges(const std::vector<unsigned int>& bulk_neighbors,
                                                              const unsigned int &component_ele_idx);
    
    friend class InspectElements;
    friend class InspectElementsAlgorithm22;
};


class InspectElementsAlgorithm22
{
public:
    InspectElementsAlgorithm22(Mesh *input_mesh);
    
    /// Runs the core algorithm for computing 2D-2D intersection in 3D.
    void compute_intersections(const std::vector<std::vector<ILpair>> &intersection_map_);
private:
    Mesh *mesh;
    
    /// Stores temporarily 2D-2D intersections.
    std::vector<IntersectionAux<2,2>> intersectionaux_storage22_;
    
    /// Object representing a triangle element A.
    Simplex<2> triaA_;
    /// Object representing a triangle element B.
    Simplex<2> triaB_;
    
    void compute_single_intersection(const ElementFullIter &eleA, const ElementFullIter &eleB);
    
    /// Auxiliary function that translates @p ElementFullIter to @p Simplex<2>.
    void update_simplex(const ElementFullIter &element, Simplex<2> & simplex);
    
    friend InspectElements;
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