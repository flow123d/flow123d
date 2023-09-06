#ifndef DG_MOCKUP_MESHES_HH_
#define DG_MOCKUP_MESHES_HH_

#include <string>
#include <vector>


/**
 * Defines and provides path to mesh files used in assembly benchmark tests.
 */

class CP {
public:
    constexpr CP(const char * str)
    : str_(str)
    {}

    const char * str_;
};

/// Define paths of test meshes here.
std::vector<std::string> meshes_table {
    "square_uniform_2_S",
	"square_uniform_2_M",
    "square_uniform_2_L"
};


#endif /* DG_MOCKUP_MESHES_HH_ */
