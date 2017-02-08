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
 * @file    inspect_elements.hh
 * @brief   Main class MixedMeshIntersections for governing intersection algorithm on meshes of combined dimensions.
 * @author  Pavel Exner
 *
 * TODO: in 2d-2d and 1d-2d check that the intersection candidate has been already computed
 *      during algorithm, not when moving to storage..
 * TODO: fix component numbering - component number is increasing weirdly...we get more components than expected
 *      then it does not work in 2d-2d
 *
 */
#ifndef INSPECT_ELEMENTS_H_
#define INSPECT_ELEMENTS_H_

#include "inspect_elements_algorithm.hh"
#include "input/input_type_forward.hh"

class Mesh; // forward declare

namespace computeintersection {

class MixedMeshIntersections;
class InspectElementsAlgorithm22;
class IntersectionLocalBase;
template<unsigned int N, unsigned int M> class IntersectionLocal;

/// Selection of intersections of different dimensions.
enum IntersectionType
{
    none = 0,
    d12 = 0x0001,
    d13 = 0x0002,
    d23 = 0x0004,
    d22 = 0x0008,
    d12_1 = 0x0010, // different algorithms for 12
    d12_2 = 0x0020,
    d12_3 = 0x0040,
    all = 0xFFFF
};

/** @brief Main class for computation of intersection of meshes of combined dimensions.
 *
 * Controls different algorithms for computation of intersection.
 * Computes intersection in the order:
 *      - 1D-3D
 *      - 2D-3D
 *      - 2D-2D (uses results from 2D-3D)
 * After computation, it transforms the internal intersection objects into proper output structure
 * for further usage.
 * The objects of different intersection types are stored in different vectors:
 * @p intersection_storageXY_ ..
 *
 * When we are on an element, we use @p element_intersections_ to get to all its intersections.
 * 
 */
class MixedMeshIntersections
{   
public:
    /// Stores 1D-3D intersections.
    std::vector<IntersectionLocal<1,3>> intersection_storage13_;
    /// Stores 2D-3D intersections.
    std::vector<IntersectionLocal<2,3>> intersection_storage23_;
    /// Stores 2D-2D intersections.
    std::vector<IntersectionLocal<2,2>> intersection_storage22_;
    /// Stores 1D-2D intersections.
    std::vector<IntersectionLocal<1,2>> intersection_storage12_;
    
    /**
     * For every element, stores list of intersections with this element.
     *
     * intersection_map_[element index][i].first = other element index
     * intersection_map_[element index][i].second = pointer to the intersection object
     */
    std::vector<std::vector<ILpair>> element_intersections_;
    
    MixedMeshIntersections(Mesh *mesh);
    ~MixedMeshIntersections();
    
    /// Calls @p InspectElementsAlgorithm<dim>, computes intersections, 
    /// move them to storage, create the map and throw away the rest.
    void compute_intersections(IntersectionType d = IntersectionType::all);
    
    // TODO: move following functions into common intersection test code.
    // Functions for tests.
    unsigned int number_of_components(unsigned int dim);
    
    /// Computes the size of the intersection in @p dim dimenstions.
    double measure_13();
    double measure_23();
    double measure_22();
    
    /// Generates a mesh file of the given name, including the intersection.
    void print_mesh_to_file_13(std::string name);
    void print_mesh_to_file_23(std::string name);
    
private:
    /// Mesh pointer.
    Mesh * mesh;
    
    InspectElementsAlgorithm<1> algorithm13_;
    InspectElementsAlgorithm<2> algorithm23_;
    InspectElementsAlgorithm22 algorithm22_;
    InspectElementsAlgorithm12 algorithm12_;
    
    /// Auxiliary function that calls InspectElementsAlgorithm<dim>.
    template<unsigned int dim> void compute_intersections(InspectElementsAlgorithm<dim> &iea,
                                                          std::vector<IntersectionLocal<dim,3>> &storage);
    void compute_intersections_22(std::vector<IntersectionLocal<2,2>> &storage);
    void compute_intersections_12(std::vector<IntersectionLocal<1,2>> &storage);
    void compute_intersections_12_2(std::vector<IntersectionLocal<1,2>> &storage);
};

    
} // END namespace


#endif // INSPECT_ELEMENTS_H_
