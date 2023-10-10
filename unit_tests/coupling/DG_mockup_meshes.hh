#ifndef DG_MOCKUP_MESHES_HH_
#define DG_MOCKUP_MESHES_HH_

#include <string>
#include <vector>


/**
 * Defines and provides path to mesh files used in assembly benchmark tests.
 */

/// Define paths of test meshes here.
std::vector<std::string> meshes_table {
    "square_2D_uniform_small",
    "square_2D_uniform_medium",
    "square_2D_uniform_big"
};


#endif /* DG_MOCKUP_MESHES_HH_ */
