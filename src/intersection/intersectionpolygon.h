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
template<unsigned int, unsigned int> class IntersectionPoint;

    /**
     * IntersectionPolygon represents a polygon from computing intersection of a triangle(Simplex<2>) and a tetrahedron(Simplex<3>)
     * Contains array of points (necessary not in order at first)
     * Flag if its "patological" (polygon is a part of side in tetrahedron)
     */
class IntersectionPolygon {

	bool pathologic_;                               ///< Pathologic flag. True if one or more ips are.
    
	std::vector<IntersectionPoint<2,3>> i_points_;  ///< Intersection points 2d-3d.
	unsigned int element_2D_idx;                    ///< Index of 2d element.
	unsigned int element_3D_idx;                    ///< Index of 3d element.

    /// Optimized polygon tracing.
    void trace_polygon_opt(std::vector<unsigned int> &prolongation_table);

    /** General polygon tracing (also resolves pathologic cases). 
     * Implements Monotone Chain algorithm (find convex hull of set of points in 2D):
     * https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
     */
    void trace_polygon_convex_hull(std::vector<unsigned int> &prolongation_table);
    
    double convex_hull_cross(const IntersectionPoint<2,3> &O,
                             const IntersectionPoint<2,3> &A,
                             const IntersectionPoint<2,3> &B) const;
                             
//     /**
//      * @return index of prolongation table for two points, which are on same line of triangle or side of tetrahedron
//      */
//     int convex_hull_prolongation_side(const IntersectionPoint<2,3> &A,
//                                       const IntersectionPoint<2,3> &B) const;
                                      
    /**
     * @return index of side of tetrahedron - if tetrahedron has polygon on it
     */
    int side_content_prolongation() const;
    
public:
	IntersectionPolygon();  ///< Default constructor.
	~IntersectionPolygon(); ///< Destructor.

    /** Constructor taking in element indices.
     */
	IntersectionPolygon(unsigned int elem2D,unsigned int elem3D);

	/** Adds new intersection point.
	 * @param i_point - intersection point to be added.
	 */
	void add_ipoint(const IntersectionPoint<2,3> &i_point);

    /// Returns true if the polygon is a pathologic case; false otherwise.
    bool is_pathologic() const;

    /// Returns intersection point of given @p index.
    const IntersectionPoint<2,3> &operator[](unsigned int index) const;

    unsigned int size() const;          ///< Returns number of intersection points.
    unsigned int ele_2d_idx() const;    ///< Returns index of 2d element.
    unsigned int ele_3d_idx() const;    ///< Returns index of 3d element.

	/** method compute local polygon area from barycentric coordinates
	 * @return double computed local area
	 * split the polygon into triangles according to the first point
	 * for every triangle compute area from barycentric coordinates
	 * 				Barycentric coords. 		Local coords.
	 * Point A		[A0,A1,A2]					[A1,A2,0]
	 * Point B	    [B0,B1,B2]					[B1,B2,0]
	 * Point C		[C0,C1,C2]					[C1,C2,0]
	 *
	 *  triangle area: (B-A)x(C-A)/2 => ((B1-A1)(C2-A2) - (B2-A2)(C1-A1))/2
	 *  => (A1(B2-C2) + B1(C2-A2) + C1(A2-B2))/2
	 *
	 *  final polygon area is sum of triangle areas
	 */
	double get_area() const;

	/**
	 * Tracing generic polygon - if flag "is_patological_" is true - traces polygon by convex hull,
	 * otherwise by optimize method
	 * it fills prolongation table
	 */
	void trace_polygon(std::vector<unsigned int> &prolongation_table);

    /// Friend output operator.
    friend std::ostream& operator<<(std::ostream& os, const IntersectionPolygon& polygon);
};

/********************************************* IMPLEMENTATION ***********************************************/

inline bool IntersectionPolygon::is_pathologic() const
{   return pathologic_; }

inline const IntersectionPoint< 2, 3 >& IntersectionPolygon::operator[](unsigned int index) const
{   ASSERT(index < i_points_.size(), "Index out of bounds.");
    return i_points_[index]; }

inline unsigned int IntersectionPolygon::size() const
{   return i_points_.size(); }

inline unsigned int IntersectionPolygon::ele_2d_idx() const
{   return element_2D_idx; }

inline unsigned int IntersectionPolygon::ele_3d_idx() const
{   return element_3D_idx; }

} // END NAMESPACE
#endif /* INTERSECTIONPOLYGON_H_ */
