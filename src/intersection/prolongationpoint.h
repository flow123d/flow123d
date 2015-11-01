#ifndef PROLONGATIONPOINT_H_
#define PROLONGATIONPOINT_H_


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

} // END NAMESPACE
#endif /* PROLONGATIONPOINT */
