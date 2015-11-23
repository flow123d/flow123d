#ifndef INTERSECTIONPOINT_H_
#define INTERSECTIONPOINT_H_

#include <armadillo>

namespace computeintersection{

//forward declare
template<unsigned int N, unsigned int M> class IntersectionPoint;
template<unsigned int N, unsigned int M> std::ostream& operator<<(std::ostream& os, const IntersectionPoint<N,M>& IP);
    
static unsigned int unset_loc_idx = -1; ///< Value for local index that is unset.
/**
 * Class represents an intersection point of simplex<N> and simplex<M>.
 * It contains barycentric coordinates of the point on both simplices.
 * Further, it contains topology information of the intersection point.
 * Namely its orientation (according to Plucker coordinates product),
 * edge and side indices, pathologic flag.
 */
template<unsigned int N, unsigned int M> class IntersectionPoint {

	arma::vec::fixed<N+1> local_coords1; ///< Barycentric coordinates of a point on simplex<N>.
	arma::vec::fixed<M+1> local_coords2; ///< Barycentric coordinates of a point on simplex<M>.

	unsigned int side_idx1; // For case N = 2, M = 3 -> index of a triangle line
	unsigned int side_idx2; // For case N = 1 or N = 2, M = 3 -> index of a tetrahedron side

	//TODO comment on what the orientation really is (relation to orientation of side,line,plucker)
	unsigned int orientation; // orientation from intersection using plucker coords

	bool is_patological_; // points is a vertex of tetrahedron or it is in side of tetrahedron or it is in edge of tetrahedron
    
	public:

    IntersectionPoint();    ///< Default constructor.
    ~IntersectionPoint(){}; ///< Destructor.
    
    /**
     * TODO can be split into two setters ? coordinate data x topolopgy data
     */
	IntersectionPoint(const arma::vec::fixed<N+1> &lc1, const arma::vec::fixed<M+1> &lc2);
	

	/// Constructor - fliping dimension of an intersectionpoint
	IntersectionPoint(IntersectionPoint<M, N> &IP);

	/* Constructor interpolates the second bary coords of IntersectionPoint<N,M-1> to IntersectionPoint<N,M>
     * Allowed only from dimension 1 to 2 and from 2 to 3
     * */
	IntersectionPoint(IntersectionPoint<N,M-1> &IP);

	/* Constructor interpolates the second bary coords of IntersectionPoint<N,M-2> to IntersectionPoint<N,M>
	 * Allowed only from dimension 1 to 3
	 * */
	IntersectionPoint(IntersectionPoint<N,M-2> &IP);

    void set_coordinates(const arma::vec::fixed<N+1> &lc1, const arma::vec::fixed<M+1> &lc2);
    void set_topology(unsigned int side1 = unset_loc_idx,
                      unsigned int side2 = unset_loc_idx,
                      unsigned int ori = 1,
                      bool patological = false);
    
//     void set_topology_SS(unsigned int edge_idx,
//                          unsigned int ori,
//                          bool pathologic);
//     void set_topology_ES(unsigned int edge_idx,
//                          unsigned int side_idx,
//                          unsigned int ori,
//                          bool pathologic);
//     void set_topology_EE(unsigned int vertex_idx,
//                          bool pathologic);

    /// Resets the object to default values.
    void clear();

    /// Returns barycentric coordinates in the Simplex<N>.
    const arma::vec::fixed<N+1> &get_local_coords1() const;
    
    /// Returns barycentric coordinates in the Simplex<M>.
    const arma::vec::fixed<M+1> &get_local_coords2() const;

    void set_side1(unsigned int s);

    void set_side2(unsigned int s);

//     void set_orientation(unsigned int o);

//     void set_is_vertex(bool iv);

//     void set_is_patological(bool ip);

    unsigned int get_side1() const;

    unsigned int get_side2() const;

    unsigned int get_orientation() const;

    bool is_vertex() const;

    bool is_patological() const;
    
	/**
	 * For convex hull polygon tracing
	 */
	bool operator<(const IntersectionPoint<N,M> &ip) const;
    
    /// Friend output operator.
    friend std::ostream& operator<< <>(std::ostream& os, const IntersectionPoint<N,M>& IP);
};


/********************************************* IMPLEMENTATION ***********************************************/

template<unsigned int N, unsigned int M>
void IntersectionPoint<N,M>::set_coordinates(const arma::vec::fixed< N + 1  >& lc1, const arma::vec::fixed< M + 1  >& lc2)
{   local_coords1 = lc1;
    local_coords2 = lc2; }

template<unsigned int N, unsigned int M>
void IntersectionPoint<N,M>::set_topology(unsigned int side1, unsigned int side2, unsigned int ori, bool pathologic)
{   side_idx1 = side1;
    side_idx2 = side2;
    orientation = ori;
    is_patological_ = pathologic;
}
    
// template<unsigned int N, unsigned int M>  
// void IntersectionPoint<N,M>::set_topology_SS(unsigned int edge_idx, unsigned int ori, bool pathologic)
// {
//     side_idx1 = unset_loc_idx;
//     side_idx2 = edge_idx;
//     orientation = ori;
//     is_patological_ = pathologic;
// }
// 
// template<unsigned int N, unsigned int M>
// void IntersectionPoint<N,M>::set_topology_ES(unsigned int edge_idx, unsigned int side_idx, unsigned int ori, bool pathologic)
// {
//     side_idx1 = edge_idx;
//     side_idx2 = side_idx;
//     orientation = ori;
//     is_patological_ = pathologic;
// }
// 
// template<unsigned int N, unsigned int M>
// void IntersectionPoint<N,M>::set_topology_EE(unsigned int vertex_idx, bool pathologic)
// {
//     side_idx1 = vertex_idx;
//     side_idx2 = unset_loc_idx;
//     orientation = 0;
//     is_patological_ = pathologic;
// }

    
template<unsigned int N, unsigned int M>
const arma::vec::fixed< N + 1  >& IntersectionPoint<N,M>::get_local_coords1() const
{   return local_coords1; }

template<unsigned int N, unsigned int M>
const arma::vec::fixed< M + 1  >& IntersectionPoint<N,M>::get_local_coords2() const
{   return local_coords2; }

template<unsigned int N, unsigned int M>
void IntersectionPoint<N,M>::set_side1(unsigned int s)
{   side_idx1 = s; }

template<unsigned int N, unsigned int M>
void IntersectionPoint<N,M>::set_side2(unsigned int s)
{   side_idx2 = s; }

// template<unsigned int N, unsigned int M>
// void IntersectionPoint<N,M>::set_orientation(unsigned int o)
// {   orientation = o; }

// template<unsigned int N, unsigned int M>
// void IntersectionPoint<N,M>::set_is_vertex(bool iv)
// {   is_vertex_ = iv; }

// template<unsigned int N, unsigned int M>
// void IntersectionPoint<N,M>::set_is_patological(bool ip)
// {   is_patological_ = ip; }

template<unsigned int N, unsigned int M>
unsigned int IntersectionPoint<N,M>::get_side1() const
{   return side_idx1; }

template<unsigned int N, unsigned int M>
unsigned int IntersectionPoint<N,M>::get_side2() const
{   return side_idx2; }

template<unsigned int N, unsigned int M>
unsigned int IntersectionPoint<N,M>::get_orientation() const
{   return orientation; }

template<unsigned int N, unsigned int M>
bool IntersectionPoint<N,M>::is_vertex() const
{   return (side_idx2 == -1); }

template<unsigned int N, unsigned int M>
bool IntersectionPoint<N,M>::is_patological() const
{   return is_patological_; }


} // END namespace
#endif /* INTERSECTIONPOINT_H_ */
