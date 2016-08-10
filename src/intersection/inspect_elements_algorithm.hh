/*!
 *
﻿* Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    intersection_local.hh
 * @brief   Classes with algorithms for computation of intersections of meshes.
 * @author  Pavel Exner
 *
 */

#ifndef INSPECT_ELEMENTS_ALGORITHM_H_
#define INSPECT_ELEMENTS_ALGORITHM_H_

#include "mesh/bounding_box.hh"
#include "mesh/mesh_types.hh"

#include "simplex.hh"

#include <queue>

class Mesh; // forward declare

namespace computeintersection {

template<unsigned int N, unsigned int M> class IntersectionPointAux;
template<unsigned int N, unsigned int M> class IntersectionAux;
template<unsigned int N, unsigned int M> class IntersectionLocal;
class IntersectionLocalBase;

class InspectElements;
class InspectElementsAlgorithm22;

/// First = element index, Second = pointer to intersection object.
typedef std::pair<unsigned int, IntersectionLocalBase*> ILpair;

/// Auxiliary function that translates @p ElementFullIter to @p Simplex<simplex_dim>.
template<unsigned int simplex_dim>
void update_simplex(const ElementFullIter &element, Simplex<simplex_dim> & simplex);
    
/** @brief Class implements algorithm for dim-dimensional intersections with 3D elements.
 * 
 * @p dim-D elements that are continuously neighbouring are called component elements.
 * 3D elements are called bulk elements.
 * 
 * Implements the initialization routine, that finds the first candidate for intersection.
 * It uses bounding boxes to fastly resolve intersection candidates. 
 * We call elements whose bounding boxes are intersecting 'candidates'.
 * 
 * Finding first candidate:
 *      - BIH search algorithm - build BIH and use it to find first candidates of components
 *      - BB search algorithm - compute only bounding boxes and search through till finding candidates
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
 * 
 * Three algorithms are implemented:
 * 
 *      - BIH only: creates BIH, uses BIH to search through all elements to find candidates
 *      - BIH search: creates BIH, uses BIH to find only first candidates of components; then uses prolongation
 *      - BB search: does not create BIH, uses bounding boxes to search through all elements to find candidates
 * 
 * Due to optimal tracing algorithm for 2d-3d, we consider tetrahedron only with positive Jacobian.
 * This is checked in assert.
 * 
 * TODO: check unit test prolongation 13d, because it has different results for BIH only and BB search
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
    
    ///@name Algorithms
    /// Runs the core algorithm for computing dimD-3D intersection.
    //@{ 
    /// Uses BIHtree to find the initial candidate of a component and then prolongates the component intersetion.
    void compute_intersections();
    /// Uses only BIHtree to find intersection candidates. (No prolongation).
    void compute_intersections_BIHtree();
    /// Tests bounding boxes intersectioss to find the initial candidate of a component
    /// and then prolongates the component intersetion. (No BIHtree).
    void compute_intersections_BB();
    //@}
    
private:
    const unsigned int undefined_elm_idx_ = -1;
    
    /// Counter for intersection among elements.
    unsigned int n_intersections_;
    unsigned int component_counter_;
    
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
    
    void assert_same_intersection(unsigned int comp_ele_idx, unsigned int bulk_ele_idx);
    
    /// A hard way to find whether the intersection of two elements has already been computed, or not.
    bool intersection_exists(unsigned int component_ele_idx, unsigned int bulk_ele_idx);
    
    /// Computes the first intersection, from which we then prolongate.
    bool compute_initial_CI(unsigned int component_ele_idx, unsigned int bulk_ele_idx, unsigned int component_idx,
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

/** @brief Implements algorithm for finding 2D-2D intersections.
 * 
 * It uses previously computed 2D-3D intersections and find candidates
 * which have intersection in the same bulk (3D) element.
 */
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
    
    /// Computes fundamental intersection of two 2D elements.
    void compute_single_intersection(const ElementFullIter &eleA, const ElementFullIter &eleB);
    
    friend InspectElements;
};


/** @brief Implements algorithm for finding 1D-2D intersections.
 * 
 * There are 3 possibilities of computation:
 * 
 * 1) 2D problem with 1D fractures
 *  2D mesh is a plane and 1D mesh is the same plane; then we can solve the intersection only in 2D,
 *  probably not using Plucker; we can apply some prolongation algorithm like in 1D-3D..
 * 
 * 2) 1D-2D intersection on 3D space
 *  1D and 2D meshes are arbitrarily placed in 3D space; then we need to compute intersections using Plucker
 *  coordinates and we cannot rely on prolongation algorithm (only if there are 2 IPs
 *  -> both lie in the same plane)
 *  - assuming that a 2D mesh component is mostly a plane in our applications (not an arbitrarily curved surface)
 * 
 * 3) 1D-2D intersection inside 3D mesh
 *  1D and 2D mesh is inside a bulk mesh (3D); if 2D-3D and 1D-3D has been computed,
 *  we can iterate over 3D intersection and get candidates if there is both 1D-3D and 2D-3D
 *  at the same time.
 *  We might need to deal with intersections outside 3D mesh, using partially algorithm (2)
 */
class InspectElementsAlgorithm12
{
public:
    InspectElementsAlgorithm12(Mesh *input_mesh);
    
    /** @brief Runs the algorithm (3): compute 1D-2D intersection inside 3D mesh.
     * It directly fills the intersection map and intersection storage.
     * The intersection map is then also used to avoid duplicit intersections.
     * @param intersection_map_ map of intersections where 1D-3D and 2D-3D must be already computed
     * @param storage vector of intersection objects 1D-2D
     */
    void compute_intersections(std::vector<std::vector<ILpair>> &intersection_map,
                               std::vector<IntersectionLocal<1,2>> &storage);
    
    /** @brief Runs the algorithm (2): compute 1D-2D intersection in 3D ambient space
     * BIH is used to find intersection candidates.
     */
    void compute_intersections_2();
private:
    Mesh *mesh;
    
    /// Stores temporarily 1D-2D intersections.
    std::vector<IntersectionAux<1,2>> intersectionaux_storage12_;
    
    /// Object representing a line element.
    Simplex<1> abscissa_;
    /// Object representing a triangle element.
    Simplex<2> triangle_;
    
    /// Computes fundamental 1D-2D intersection of candidate pair.
//     void compute_single_intersection(const ElementFullIter &comp_ele, const ElementFullIter &bulk_ele);
    
    friend InspectElements;
};

} // END namespace

#endif // INSPECT_ELEMENTS_ALGORITHM_H_