#ifndef DG_MOCKUP_MESHES_HH_
#define DG_MOCKUP_MESHES_HH_

#include <string>
#include <vector>


/**
 * Defines and provides path to mesh files used in assembly benchmark tests.
 */

/// Define paths of test meshes here.
std::vector<std::string> meshes_table {
    "square_2D_uniform",
    "square_2D_refined",
    "lshape_2D_uniform",
    "lshape_2D_refined",
    "cube_3D_uniform",
    "cube_3D_refined",
    "lshape_3D_uniform",
    "lshape_3D_refined"
};

/// Define sizes of meshes here.
std::vector<std::string> meshes_sizes {
    "small",
    "medium",
    "big"
};


#endif /* DG_MOCKUP_MESHES_HH_ */
