#ifndef DG_MOCKUP_MESHES_HH_
#define DG_MOCKUP_MESHES_HH_

#include <string>
#include <vector>


/**
 * Defines and provides path to mesh files used in assembly benchmark tests.
 */

/// Define paths of test meshes here.
std::vector<std::string> meshes_table {
    "square_2D_small_uniform",
    "square_2D_medium_uniform",
    "square_2D_big_uniform",
    "square_2D_small_refined",
    "square_2D_medium_refined",
    "square_2D_big_refined",
    "lshape_2D_small_uniform",
    "lshape_2D_medium_uniform",
    "lshape_2D_big_uniform",
    "lshape_2D_small_refined",
    "lshape_2D_medium_refined",
    "lshape_2D_big_refined"/*,
    "cube_3D_small_uniform",
    "cube_3D_medium_uniform",
    "cube_3D_big_uniform",
    "cube_3D_small_refined",
    "cube_3D_medium_refined",
    "cube_3D_big_refined",
    "lshape_3D_small_uniform",
    "lshape_3D_medium_uniform",
    "lshape_3D_big_uniform",
    "lshape_3D_small_refined",
    "lshape_3D_medium_refined",
    "lshape_3D_big_refined"*/
};


#endif /* DG_MOCKUP_MESHES_HH_ */
