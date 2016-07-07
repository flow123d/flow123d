/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
 * @file    tracing_algorithm.hh
 * @brief   Implementation of tracing algorithms for polygonal intersections.
 * @author  Viktor Fris, Pavel Exner
 *
 */

#ifndef TRACING_ALGORITHM_H_
#define TRACING_ALGORITHM_H_

#include "system/global_defs.h"

namespace computeintersection{

//forwward declare
template<unsigned int, unsigned int> class IntersectionPointAux;
template<unsigned int, unsigned int> class IntersectionAux;

/** @brief Static class implementing tracing algorithms for polynomial intersection.
 */
class Tracing
{
public:
    /** @brief Trace general polygon.
     * Decides which algorithm will be used, according to type of polygon: 
     *      - pathologic -> use general convex hull algorithm
     *      - not pathologic -> use optimized algorithm exploiting Plucker data
     * It fills the prolongation table for later usage in prolongation routines.
     */
    static void trace_polygon(std::vector<unsigned int> &prolongation_table, IntersectionAux<2,3> &polygon);

private:
    /// Optimized polygon tracing.
    static void trace_polygon_opt(std::vector<unsigned int> &prolongation_table, IntersectionAux<2,3> &polygon);

    /** General polygon tracing (also resolves pathologic cases). 
     * Implements Monotone Chain algorithm (find convex hull of set of points in 2D):
     * https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
     */
    static void trace_polygon_convex_hull(std::vector<unsigned int> &prolongation_table, IntersectionAux<2,3> &polygon);
    
    /** 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
     * @return positive value, if OAB makes a counter-clockwise turn;
     *      negative for clockwise turn, and zero if the points are collinear.
     */
    static double convex_hull_cross(const IntersectionPointAux<2,3> &O,
                                    const IntersectionPointAux<2,3> &A,
                                    const IntersectionPointAux<2,3> &B);
                             
    /** @brief Auxiliary function for prolongation table construction
     * Counts the IPs on each face of tetrahedron.
     * If there is a side with more than 2 IPs on it, it means, that triangle and face are in the same plane.
     * Then index of this side is returned and marks that we do not have to prolongate over this face.
     * @return tetrahedron face index or -1
     */
    static int side_content_prolongation(IntersectionAux<2,3> &polygon);
};


} // END NAMESPACE
#endif // TRACING_ALGORITHM_H_
