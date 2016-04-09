/*
 * IntersectionPolygon.h
 *
 *  Created on: 8.4.2015
 *      Author: viktor
 */

#ifndef INTERSECTIONPOLYGON_H_
#define INTERSECTIONPOLYGON_H_

#include "system/system.hh"

namespace computeintersection{

//forwward declare
template<unsigned int, unsigned int> class IntersectionPointAux;
template<unsigned int, unsigned int> class IntersectionAux;

class Tracing
{
public:
    /**
     * Tracing generic polygon - if flag "is_patological_" is true - traces polygon by convex hull,
     * otherwise by optimize method
     * it fills prolongation table
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
    
    static double convex_hull_cross(const IntersectionPointAux<2,3> &O,
                                    const IntersectionPointAux<2,3> &A,
                                    const IntersectionPointAux<2,3> &B);
                             
    /**
     * @return index of side of tetrahedron - if tetrahedron has polygon on it
     */
    static int side_content_prolongation(IntersectionAux<2,3> &polygon);
};


} // END NAMESPACE
#endif /* INTERSECTIONPOLYGON_H_ */
