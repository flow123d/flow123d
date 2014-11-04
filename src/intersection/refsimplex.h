#include <armadillo>
#include <array>
namespace computeintersection{
/*
 * Ordering of nodes and sides in reference simplices
 * =================================================
 *
 * 1D simplex (line segment)   2D simplex (triangle)        3D simplex (tetrahedron)
 *
 *                                                                            y
 *                                                                          .
 *                                                                        ,/
 *                                                                       /
 *                                                                    3
 *                             y                                    ,/|`\
 *                             ^                                  ,/  |  `\
 *                             |                                ,/    '.   `\
 *                             2                              ,/       |     `\
 *                             |`\                          ,/         |       `\
 *                             |  `\                       0-----------'.--------1 --> x
 *                             |    `\                      `\.         |      ,/
 *                             |      `\                       `\.      |    ,/
 *                             |        `\                        `\.   '. ,/
 * 0----------1 --> x          0----------1 --> x                    `\. |/
 *                                                                      `2
 *                                                                         `\.
 *                                                                            ` z
 *
 * side id  node ids           side id  node ids           side id  node ids   normal
 * 0        0                  0        0,1                0        0,1,2      OUT
 * 1        1                  1        0,2                1        0,1,3      IN
 *                             2        1,2                2        0,2,3      OUT
 *                                                         3        1,2,3      IN
 *
 *
 */

template<unsigned int dim>
class RefSimplex
{
public:

	template<unsigned int subdim> static RefSimplex<subdim> SubRefSimplex();


	template<unsigned int subdim> inline static std::array<arma::vec::fixed<dim+1>,subdim+1> bary_coords(unsigned int sid){
		    //ASSERT(subdim < dim, "Sub-dimension is bigger than dimension!");
			//xprintf(Msg, "barycoods \n");
			std::array<arma::vec::fixed<dim+1>,subdim+1> bary_c;


			for(unsigned int i = 0; i < subdim+1; i++){
				if((dim-subdim) == 2){
					bary_c[i] = RefSimplex<dim>::node_coords(RefSimplex<dim>::line_nodes[sid][i]);
				}else{
					bary_c[i] = RefSimplex<dim>::node_coords(RefSimplex<dim>::side_nodes[sid][i]);
				}
				//bary_c[i] = RefSimplex<dim>::node_coords(RefSimplex<dim>::line_nodes[sid][i] : RefSimplex<dim>::side_nodes[sid][i]);
			}

			return bary_c;
	};

	template<unsigned int subdim> inline static arma::vec::fixed<dim+1> interpolate(arma::vec::fixed<subdim+1> coord, int sub_simplex_idx){

		std::array<arma::vec::fixed<dim+1>, subdim+1> simplex_M_vertices = RefSimplex<dim>::bary_coords<subdim>(sub_simplex_idx);
		arma::vec::fixed<dim+1> sum;
		sum.zeros();
		for(int i=0; i<subdim+1; i++) sum += coord[i]*simplex_M_vertices[i];
		return sum;
	};

	inline static arma::vec::fixed<dim+1> line_barycentric_interpolation(arma::vec::fixed<dim+1> first_coords, arma::vec::fixed<dim+1> second_coords, double first_theta, double second_theta, double theta){

		//cout << "Barycentric interpolation (first theta, theta, second theta) - (" << first_theta << "," << theta << "," << second_theta << ")" << endl;

		arma::vec::fixed<dim+1> bary_interpolated_coords;

		bary_interpolated_coords = ((theta - first_theta) * second_coords + (second_theta - theta) * first_coords)/(second_theta - first_theta);

		return bary_interpolated_coords;
	};

	//inline static arma::vec::fixed<dim+1> point_interpolation(arma::vec3 &point_coords, Simplex<3> &tetrahedron){};

	/**
	 * Return barycentric coordinates of given node.
	 * @param nid Node number.
	 */
	static arma::vec::fixed<dim+1> node_coords(unsigned int nid);

	/**
	 * Compute normal vector to a given side.
	 * @param sid Side number.
	 */
	static arma::vec::fixed<dim> normal_vector(unsigned int sid);

	/**
	 * Return index of 1D line, shared by two faces @p f1 and @p f2 of the reference tetrahedron.
	 * Implemented only for @p dim == 3.
	 */
	static unsigned int line_between_faces(unsigned int f1, unsigned int f2);

	/**
	 * Number of sides.
	 */
	static const unsigned int n_sides = dim + 1;

	/**
	 * Number of vertices.
	 */
	static const unsigned int n_nodes = dim + 1;

	/**
	 * Number of nodes on one side.
	 */
	static const unsigned int n_nodes_per_side = dim;

	/// Number of lines on boundary of one side.
	static const unsigned int n_lines_per_side = (unsigned int)((dim * (dim - 1)) / 2);//( dim == 3 ? 3 : 0);// Kombinační číslo dim nad dvěma

	/// Number of lines, i.e. @p object of dimension @p dim-2 on the boundary of the reference element.
	static const unsigned int n_lines = (unsigned int)((dim * (dim + 1)) / 2); //( dim == 3 ? 6 : dim == 2 ? 3 : dim == 1 ? 1 : 0); součet posloupnosti

	/**
	 * Node numbers for each side.
	 */
	static const unsigned int side_nodes[n_sides][n_nodes_per_side];

	/**
	 * Indices of 1D lines of the 2D sides of an tetrahedron. Nonempty only for @p dim==3.
	 */
	static const unsigned int side_lines[n_sides][n_lines_per_side];

	/**
	 * Nodes of 1D lines of the tetrahedron.
	 */
    static const unsigned int line_nodes[n_lines][2];



	static const unsigned int line_sides[n_lines][2];


    /**
	 * Number of permutations of nodes on sides.
	 * dim   value
	 * -----------
	 * 1     1
	 * 2     2
	 * 3     6
	 */
	static const unsigned int n_side_permutations = (dim+1)*(2*dim*dim-5*dim+6)/6;

	/**
	 * Permutations of nodes on sides.
	 */
	static const unsigned int side_permutations[n_side_permutations][n_nodes_per_side];

	/**
	 * For a given permutation @p p of nodes finds its index within @p side_permutations.
	 * @param p Permutation of nodes.
	 */
	static unsigned int permutation_index(unsigned int p[n_nodes_per_side]);

};
} // END namespace
