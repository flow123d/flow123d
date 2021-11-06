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

// #include "simplex.hh"

#include <queue>

class Mesh; // forward declare
class BIHTree;

template<unsigned int N, unsigned int M> class IntersectionPointAux;
template<unsigned int N, unsigned int M> class IntersectionAux;
template<unsigned int N, unsigned int M> class IntersectionLocal;
class IntersectionLocalBase;
class InspectElementsAlgorithm22;
template <int spacedim> class ElementAccessor;


/// First = element index, Second = pointer to intersection object.
typedef std::pair<unsigned int, IntersectionLocalBase*> ILpair;

template<unsigned int dimA, unsigned int dimB>
class IntersectionAlgorithmBase{
public:
    IntersectionAlgorithmBase(Mesh* mesh);
protected:
   
    /// Auxiliary function that translates @p ElementAccessor<3> to @p Simplex<simplex_dim>.
//    template<unsigned int simplex_dim>
//    void update_simplex(const ElementAccessor<3> &element, Simplex<simplex_dim> & simplex);
    
    /// Mesh pointer.
    Mesh *mesh;
    
    const unsigned int undefined_elm_idx_ = -1;
    
    /// Objects representing single elements.
//     Simplex<dimA> simplexA;
//     Simplex<dimB> simplexB;
};

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
class InspectElementsAlgorithm : public IntersectionAlgorithmBase<dim,3>
{
public:  
    InspectElementsAlgorithm(Mesh *_mesh);
    ~InspectElementsAlgorithm();
    
    ///@name Algorithms
    /// Runs the core algorithm for computing dimD-3D intersection.
    //@{ 
    /// Uses BIHtree to find the initial candidate of a component and then prolongates the component intersetion.
    void compute_intersections(const BIHTree& bih);
    /// Uses only BIHtree to find intersection candidates. (No prolongation).
    void compute_intersections_BIHtree(const BIHTree& bih);
    /// Tests bounding boxes intersectioss to find the initial candidate of a component
    /// and then prolongates the component intersetion. (No BIHtree).
    void compute_intersections_BB();
    //@}
    
private:
    using IntersectionAlgorithmBase<dim,3>::mesh;
    using IntersectionAlgorithmBase<dim,3>::undefined_elm_idx_;
//     using IntersectionAlgorithmBase<dim,3>::simplexA;
//     using IntersectionAlgorithmBase<dim,3>::simplexB;
    
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
    
    /// Counter for intersection among elements.
    unsigned int n_intersections_;
    
    /// Prolongation queue in the component mesh.
    std::queue<Prolongation> component_queue_;
    /// Prolongation queue in the bulk mesh.
    std::queue<Prolongation> bulk_queue_;
    
    // Array of flags, which elements are computed
    std::vector<bool> closed_elements;
    std::vector<unsigned int> last_slave_for_3D_elements;
    
    /// Elements bounding boxes.
    std::vector<BoundingBox> elements_bb;
    /// Bounding box of all 3D elements.
    BoundingBox mesh_3D_bb;
    
    /// Resulting vector of intersections.
    std::vector<std::vector<IntersectionAux<dim,3>>> intersection_list_;
    
    /// Initialization.
    /// Sets vector sizes and computes bulk bounding box.
    void init();
    
    /// Computes bounding boxes of all elements. Fills @p elements_bb and @p mesh_3D_bb.
    void compute_bounding_boxes();
    
    void assert_same_intersection(unsigned int comp_ele_idx, unsigned int bulk_ele_idx);
    
    /// A hard way to find whether the intersection of two elements has already been computed, or not.
    bool intersection_exists(unsigned int component_ele_idx, unsigned int bulk_ele_idx);
    
    /// Computes the first intersection, from which we then prolongate.
    bool compute_initial_CI(const ElementAccessor<3> &comp_ele, const ElementAccessor<3> &bulk_ele);
    
    /// Finds neighbouring elements that are new candidates for intersection and pushes
    /// them into component queue or bulk queue.
    void prolongation_decide(const ElementAccessor<3> &comp_ele, const ElementAccessor<3> &bulk_ele,
                             IntersectionAux<dim,3> is);
    
    /// Computes the intersection for a candidate in a queue and calls @p prolongation_decide again.
    void prolongate(const Prolongation &pr);

    template<unsigned int ele_dim>
    std::vector< unsigned int > get_element_neighbors(const ElementAccessor<3>& ele,
                                                      unsigned int ip_dim,
                                                      unsigned int ip_obj_idx);
    
//     unsigned int create_prolongation_to_bulk(unsigned int bulk_ele_idx,
//                                              unsigned int component_ele_idx);
    
    unsigned int create_prolongation(unsigned int bulk_ele_idx,
                                     unsigned int component_ele_idx,
                                     std::queue<Prolongation>& queue);
    
    friend class MixedMeshIntersections;
};

/** @brief Implements algorithm for finding 2D-2D intersections.
 * 
 * It uses previously computed 2D-3D intersections and find candidates
 * which have intersection in the same bulk (3D) element.
 */
class InspectElementsAlgorithm22 : public IntersectionAlgorithmBase<2,2>
{
public:
    InspectElementsAlgorithm22(Mesh *input_mesh);
    
    /// Runs the core algorithm for computing 2D-2D intersection in 3D.
//     void compute_intersections(const std::vector<std::vector<ILpair>> &intersection_map_);
    
    void compute_intersections(std::vector<std::vector<ILpair>> &intersection_map,
                               std::vector<IntersectionLocal<2,2>> &storage);
    
private:
    unsigned int component_counter_;
    const unsigned int unset_comp = (unsigned int)(-1);
    
    /// Stores temporarily 2D-2D intersections.
    std::vector<IntersectionAux<2,2>> intersectionaux_storage22_;
    
    // Auxiliary vector which keeps component indices for 2d elements.
    std::vector<unsigned int> component_idx_;
    
    /// Computes fundamental intersection of two 2D elements.
    void compute_single_intersection(const ElementAccessor<3> &eleA, const ElementAccessor<3> &eleB,
                                     std::vector<IntersectionLocal<2,2>> &storage);
    
    /// Creates numbering of the 2D components and fills component_idx_ vector.
    void create_component_numbering();
    
    /// Auxiliary function for front-advancing alg. for component numbering.
    //void prolongate22(const ElementAccessor<3>& ele, std::queue<unsigned int>& queue);
    
    friend class MixedMeshIntersections;
};


/** @brief Implements algorithm for finding 1D-2D intersections.
 * 
 * There are 3 possibilities of computation:
 * 
 * 1) 2D problem with 1D fractures
 *  2D mesh is a plane and 1D mesh is the same plane; then we can solve the intersection only in 2D,
 *  probably not using Plucker; we can apply some prolongation algorithm like in 1D-3D..
 * 
 * 2) 1D-2D intersection in 3D space
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
class InspectElementsAlgorithm12 : public IntersectionAlgorithmBase<1,2>
{
public:
    InspectElementsAlgorithm12(Mesh *input_mesh);
    
    /** @brief Runs the algorithm (3): compute 1D-2D intersection inside 3D mesh.
     * It directly fills the intersection map and intersection storage.
     * The intersection map is then also used to avoid duplicit intersections.
     * @param intersection_map_ map of intersections where 1D-3D and 2D-3D must be already computed
     * @param storage vector of intersection objects 1D-2D
     */
    void compute_intersections_3(std::vector<std::vector<ILpair>> &intersection_map,
                               std::vector<IntersectionLocal<1,2>> &storage);
    
    /** @brief Runs the algorithm (2): compute 1D-2D intersection in 3D ambient space
     * BIH is used to find intersection candidates.
     */
    void compute_intersections_2(const BIHTree& bih);
    
    /** @brief Runs the algorithm (1): compute 1D-2D intersection in 2D plane.
     * BIH is used to find intersection candidates.
     */
    void compute_intersections_1(const BIHTree& bih);
    
private: 
    /// Stores temporarily 1D-2D intersections.
    std::vector<IntersectionAux<1,2>> intersectionaux_storage12_;
    
    /// Computes fundamental 1D-2D intersection of candidate pair.
//     void compute_single_intersection(const ElementAccessor<3> &comp_ele, const ElementAccessor<3> &bulk_ele);
    
    friend class MixedMeshIntersections;
};


#endif // INSPECT_ELEMENTS_ALGORITHM_H_
