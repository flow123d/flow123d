#ifndef PROLONGATION_H_
#define PROLONGATION_H_


namespace computeintersection{

/**
 * Simple class defines indices of elements for later processing of computing intersection 1D-3D
 * TODO:
 * - is the old index needed?
 */
struct ProlongationPoint{
    unsigned int elm_1D_idx;        ///< Index of 1d element.
    unsigned int elm_3D_idx;        ///< Index of 3d element.
    unsigned int elm_3D_idx_old;    ///< Index of the old 3d element.
};

/**
 * Simple class defines indices of elements for later processing of computing intersection 2D-3D
 * 
 * TODO:
 * - comment
 * - is the old index needed? if not, then the comment below in redundant
 * - split prolongationline into Prolongation2D (neighbor is 2D) and Prolongation3D (neighbor is 3D)
 */

struct ProlongationLine{
    unsigned int elm_2D_idx;
    unsigned int elm_3D_idx;
    unsigned int dictionary_idx; // index to dictionary with all intersections associated with index of 2D element
    int elm_2D_idx_old;
    int elm_3D_idx_old;
};

} // END NAMESPACE
#endif /* PROLONGATION */
