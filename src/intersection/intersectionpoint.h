#ifndef INTERSECTIONPOINT_H_
#define INTERSECTIONPOINT_H_

#include <armadillo>

namespace computeintersection{

//forward declare
template<unsigned int N, unsigned int M> class IntersectionPoint;
template<unsigned int N, unsigned int M> std::ostream& operator<<(std::ostream& os, const IntersectionPoint<N,M>& IP);
    
static const unsigned int unset_loc_idx = -1; ///< Value for local index that is unset.
/**
 * Class represents an intersection point of simplex<N> and simplex<M>.
 * It contains barycentric coordinates of the point on both simplices.
 * Further, it contains topology information of the intersection point.
 * Namely its orientation (according to Plucker coordinates product),
 * edge and side indices, pathologic flag.
 */
template<unsigned int N, unsigned int M> class IntersectionPoint {
    
	arma::vec::fixed<N+1> local_coords1_; ///< Barycentric coordinates of a point on simplex<N>.
	arma::vec::fixed<M+1> local_coords2_; ///< Barycentric coordinates of a point on simplex<M>.

	unsigned int side_idx1_; // For case N = 2, M = 3 -> index of a triangle line
	unsigned int side_idx2_; // For case N = 1 or N = 2, M = 3 -> index of a tetrahedron side

	//TODO comment on what the orientation really is (relation to orientation of side,line,plucker)
	unsigned int orientation_; // orientation from intersection using plucker coords

	bool pathologic_; // points is a vertex of tetrahedron or it is in side of tetrahedron or it is in edge of tetrahedron
	
	unsigned int dim_A_;    ///< Dimension of the object A of intersection.
    unsigned int dim_B_;    ///< Dimension of the object B of intersection.
    
public:

    IntersectionPoint();    ///< Default constructor.
    ~IntersectionPoint(){}; ///< Destructor.
    
    /**
     * Constructor taking barycentric coordinates on simplices as input parameters.
     * @param lc1 - barycentric coordinates in Simplex<N>
     * @param lc2 - barycentric coordinates in Simplex<M>
     */
	IntersectionPoint(const arma::vec::fixed<N+1> &lc1, const arma::vec::fixed<M+1> &lc2,
                      unsigned int dim_A = N, unsigned int dim_B = M);
	

	/// Constructor - fliping dimension of an intersection point.
	IntersectionPoint(IntersectionPoint<M, N> &IP);

	/** Constructor interpolates the second bary coords of IntersectionPoint<N,M-1> to IntersectionPoint<N,M>
     * Allowed only from dimension 1 to 2 and from 2 to 3.
     * @param  side_idx2 is the index of object 2 of IntersectionPoint<N,M-1> in object 2 of IntersectionPoint<N,M>
     * */
	IntersectionPoint(IntersectionPoint<N,M-1> &IP, unsigned int side_idx2);

	/** Constructor interpolates the second bary coords of IntersectionPoint<N,M-2> to IntersectionPoint<N,M>
	 * Allowed only from dimension 1 to 3.
     * @param  side_idx2 is the index of object 2 of IntersectionPoint<N,M-2> in object 2 of IntersectionPoint<N,M>
	 * */
	IntersectionPoint(IntersectionPoint<N,M-2> &IP, unsigned int side_idx2);

    /// Setter for coordinates.
    void set_coordinates(const arma::vec::fixed<N+1> &lc1, const arma::vec::fixed<M+1> &lc2);
    
    /// Setter for topology data.
    void set_topology(unsigned int side1,// = unset_loc_idx,
                      unsigned int dim_A,// = N,
                      unsigned int side2,// = unset_loc_idx,
                      unsigned int dim_B// = M,
                     );
    
    /// Setter for Plucker flags.
    void set_plucker_flags(unsigned int ori, bool pathologic);
    
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
    const arma::vec::fixed<N+1> &local_coords1() const;
    
    /// Returns barycentric coordinates in the Simplex<M>.
    const arma::vec::fixed<M+1> &local_coords2() const;

    void set_side1(unsigned int idx);  ///<  Sets the index of Simplex<N>.
    void set_side2(unsigned int idx);   ///<  Sets the index of Simplex<M>.

//     void set_orientation(unsigned int o);

//     void set_is_vertex(bool iv);

//     void set_is_patological(bool ip);

    unsigned int dim_A() const;         ///< Returns dimension of object A.
    unsigned int dim_B() const;         ///< Returns dimension of object B.
    unsigned int side_idx1() const;     ///<  Returns the index of Simplex<N>.
    unsigned int side_idx2() const;     ///<  Returns the index of Simplex<M>.
    unsigned int orientation() const;   ///<  Returns the orientation.
    bool is_pathologic() const;         ///<  Returns true, if this is a pathologic case.
    
    /**
     * Returns true if the IP is a vertex (of a triangle).
     * Is determined according to topology data:
     *      side1        side2        ori
     * S-S    -         edge idx       *
     * E-S  edge idx    side idx       *
     * E-E  vertex idx     -           *
     * 
     * where '-' means unset, '*' means a value
     */
    bool is_vertex() const;
    
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
{   local_coords1_ = lc1;
    local_coords2_ = lc2; }

template<unsigned int N, unsigned int M>
void IntersectionPoint<N,M>::set_topology(unsigned int side1,unsigned int dim_A, unsigned int side2, unsigned int dim_B)
{   side_idx1_ = side1;
    side_idx2_ = side2;
    dim_A_ = dim_A;
    dim_B_ = dim_B;
}
    
template<unsigned int N, unsigned int M>
void IntersectionPoint<N,M>::set_plucker_flags(unsigned int ori, bool pathologic)
{
    orientation_ = ori;
    pathologic_ = pathologic;
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
unsigned int IntersectionPoint<N,M>::dim_A() const
{   return dim_A_; }

template<unsigned int N, unsigned int M>
unsigned int IntersectionPoint<N,M>::dim_B() const
{   return dim_B_; }

template<unsigned int N, unsigned int M>
const arma::vec::fixed< N + 1  >& IntersectionPoint<N,M>::local_coords1() const
{   return local_coords1_; }

template<unsigned int N, unsigned int M>
const arma::vec::fixed< M + 1  >& IntersectionPoint<N,M>::local_coords2() const
{   return local_coords2_; }

template<unsigned int N, unsigned int M>
void IntersectionPoint<N,M>::set_side1(unsigned int idx)
{   side_idx1_ = idx; }

template<unsigned int N, unsigned int M>
void IntersectionPoint<N,M>::set_side2(unsigned int idx)
{   side_idx2_ = idx; }

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
unsigned int IntersectionPoint<N,M>::side_idx1() const
{   return side_idx1_; }

template<unsigned int N, unsigned int M>
unsigned int IntersectionPoint<N,M>::side_idx2() const
{   return side_idx2_; }

template<unsigned int N, unsigned int M>
unsigned int IntersectionPoint<N,M>::orientation() const
{   return orientation_; }

template<unsigned int N, unsigned int M>
bool IntersectionPoint<N,M>::is_vertex() const
{   return (side_idx2_ == -1); }

template<unsigned int N, unsigned int M>
bool IntersectionPoint<N,M>::is_pathologic() const
{   return pathologic_; }


} // END namespace
#endif /* INTERSECTIONPOINT_H_ */
