    /*
 * observe.hh
 *
 *  Created on: Jun 28, 2016
 *      Author: jb
 */

#ifndef SRC_IO_OBSERVE_HH_
#define SRC_IO_OBSERVE_HH_

#include <boost/exception/info.hpp>          // for error_info::~error_info<...
#include <iosfwd>                            // for ofstream, ostream
#include <map>                               // for map, map<>::value_compare
#include <memory>                            // for shared_ptr
#include <new>                               // for operator new[]
#include <string>                            // for string, operator<<
#include <vector>                            // for vector
#include <armadillo>
#include "input/accessors.hh"                // for Array (ptr only), Record
#include "input/input_exception.hh"          // for DECLARE_INPUT_EXCEPTION
#include "system/exceptions.hh"              // for operator<<, ExcStream, EI
#include "system/armadillo_tools.hh"         // for Armadillo vec string
#include "system/index_types.hh"             // for LongIdx
#include "mesh/range_wrapper.hh"
#include "tools/general_iterator.hh"
#include "la/distribution.hh"

class ElementDataCacheBase;
class Mesh;
namespace Input { namespace Type { class Record; } }
template <typename T> class ElementDataCache;



/**
 * Helper class stores base data of ObservePoint and allows to evaluate
 * the nearest point to input_point_.
 */
class ObservePointData {
public:
	/// Constructor
	ObservePointData()
	: distance_(numeric_limits<double>::infinity()) {};

	/// Final element of the observe point. The index in the mesh.
	unsigned int element_idx_;

	/// Global coordinates of the observation point.
	arma::vec3 global_coords_;

	/// Local (barycentric) coordinates on the element.
	arma::vec local_coords_;

	/// Distance of found projection from the initial point.
	/// If we find more candidates we pass in the closest one.
	double distance_;

	/// Actual process of the observe point.
	unsigned int proc_;

	/// Global index of the observe point.
	LongIdx global_idx_;

	/// Local index on actual process of the observe point.
	LongIdx local_idx_;
};


/**
 * Class representing single observe point, used internally by the class Observe.
 * Members: input_pos_, snap_dim_, snap_region_name_ are set in constructor. Should be checked before passed in.
 * Members: element_idx_, global_coords_, local_coords_ are derived, set in Observe::find_observe_points.
 */
class ObservePoint {
public:
    DECLARE_INPUT_EXCEPTION(ExcNoInitialPoint,
            << "Failed to find the element containing the initial observe point.\n");
    TYPEDEF_ERR_INFO(EI_RegionName, std::string);
    TYPEDEF_ERR_INFO(EI_PointName, std::string);
    TYPEDEF_ERR_INFO(EI_Point, arma::vec3);
    TYPEDEF_ERR_INFO(EI_ClosestEle, ObservePointData);
    DECLARE_INPUT_EXCEPTION(ExcNoObserveElementCandidates,
            << "Failed to find any element in the search radius of the observe point " << EI_PointName::qval
            << " with given coordinates " << field_value_to_yaml(EI_Point::ref(*this)) << ".\n"
            << "The closest element has index " << EI_ClosestEle::ref(*this).element_idx_ << ", its distance is " << EI_ClosestEle::ref(*this).distance_ << ".\n"
            << "Solution: check the position of the observe point, possibly increase the maximal snapping distance "
            << "(keys: observe_points:search_radius, mesh:global_snap_radius)"<< "\n");
    DECLARE_INPUT_EXCEPTION(ExcNoObserveElement,
            << "Failed to find any element in the search radius of the observe point"  << EI_PointName::qval
            << " inside the snap region: " << EI_RegionName::qval << ".\n"
            << "The observe point coordinates are " << field_value_to_yaml(EI_Point::ref(*this)) << ".\n"
            << "The closest element (outside the snap region) has index " << EI_ClosestEle::ref(*this).element_idx_
            << ", its distance is " << EI_ClosestEle::ref(*this).distance_ << ".\n"
            << "Solution: check the position/region of the observe point, possibly increase the maximal snapping distance "
            << "(keys: observe_points:search_radius, mesh:global_snap_radius)"<< "\n");

    static const Input::Type::Record & get_input_type();

    /**
     * Return index of observation point in the mesh.
     */
    inline unsigned int element_idx() const
    { return observe_data_.element_idx_; }

    /**
     * Return global coordinates of the observation point.
     */
    inline arma::vec3 global_coords() const
    { return observe_data_.global_coords_; }

protected:
    /**
     *  Default constructor just for testing.
     */
    ObservePoint();

    /**
     * Constructor. Read from input.
     */
    ObservePoint(Input::Record in_rec, Mesh &mesh, unsigned int point_idx);

    /**
     * Returns true if we have already found any observe element.
     */
    bool have_observe_element();

    /**
     *  Snap to the center of closest subelement with dimension snap_dim_.
     *  This makes final adjustment of global_coords_ and local_coords_.
     */
    void snap(Mesh &mesh);


    /**
     * Find the observe element and the definitive observe point.
     *
     * Algorithm:
     * 1. find element containing the point (initial element)
     * 2. check initial element for region match possibly set it as (observe element)
     * 3. add neighbours into processing_list for the next level
     * 4. while we have no observe element: pass through the processing list of the current level
     * 5.   if element match the region, project and clip the init point, update observe element.
     * 6. snapping on the observe element
     *
     */
    void find_observe_point(Mesh &mesh);

    /**
     * Output the observe point information into a YAML formated stream, indent by
     * given number of spaces + "- ".
     */
    void output(ostream &out, unsigned int indent_spaces, unsigned int precision);

    /// Project point to given element by dimension of this element.
    ObservePointData point_projection(unsigned int i_elm, ElementAccessor<3> elm);

    /// Index in the input array.
    Input::Record in_rec_;

    /// Observation point name.
    std::string name_;

	/**
	 * Snap to the center of the object of given dimension.
	 * Value 4 and greater means no snapping.
	 */
	unsigned int snap_dim_;

    /**
     * Region of the snapping element.
     */
    string snap_region_name_;

    /**
     * Maximal distance of observe element from input point.
     */
    double max_search_radius_;

	/// Input coordinates of the initial position of the observation point.
	arma::vec3 input_point_;

    /// Helper object stored projection data
    ObservePointData observe_data_;

	/// Only Observe should use this class directly.
    friend class Observe;

};


class ObservePointAccessor;

/**
 * This class takes care about the observe points in the output stream, storing observe values of the fields and
 * their output in the YAML format.
 */
class Observe {
public:

    typedef std::shared_ptr<ElementDataCacheBase> OutputDataPtr;
    typedef std::map< string,  OutputDataPtr > OutputDataFieldMap;

    /**
     * Construct the observation object.
     *
     * observe_name - base name of the output file, the equation name.
     * mesh - the mesh used for search for the observe points
     * in_array - the array of observe points
     */
    Observe(string observe_name, Mesh &mesh, Input::Array in_array, unsigned int precision, std::string unit_str);

    /// Destructor, must close the file.
    ~Observe();


    /**
     * Provides a vector of element indices on which the observation values are evaluated.
     * This can be used to evaluate derived fields only on these elements in the times not selected to
     * full output.
     */
    inline const std::vector<unsigned int> &observed_elements() const
            { return observed_element_indices_;}

    /**
     * Output file header.
     */
    void output_header();

    /**
     * Sets next output time frame of observe. If the table is full, writes field values to the output
     * file. Using the YAML format. Argument flush starts writing to output file explicitly.
     */
    void output_time_frame(bool flush);

    /**
     * Return \p points_ vector
     */
    inline const std::vector<ObservePoint> & points() const
    { return points_; }

    /**
     * Return point distribution
     */
    inline const Distribution * point_ds() const
    { return point_ds_; }

    /// Returns local range of observe points
    Range<ObservePointAccessor> local_range() const;

    /**
     * Prepare data for computing observe values.
     *
     * Method:
     *  - check that all fields of one time frame are evaluated at the same time
     *  - find and return ElementDataCache of given field_name, create its if doesn't exist.
     *
     * @param field_name Quantity name of founding ElementDataCache
     * @param field_time Actual computing time
     * @param n_rows     Count of rows of data cache (used only if new cache is created)
     * @param n_cols     Count of columns of data cache (used only if new cache is created)
     */
    template <typename T>
    ElementDataCache<T> & prepare_compute_data(std::string field_name, double field_time, unsigned int n_rows, unsigned int n_cols);



protected:
    /// Maximal size of observe values times vector
    static const unsigned int max_observe_value_time;

    // MPI rank.
    int rank_;

    /// Full information about observe points.
    std::vector<ObservePoint> points_;
    /// Elements of the o_points.
    std::vector<unsigned int> observed_element_indices_;

    /// Stored field values.
    OutputDataFieldMap observe_field_values_;


    /// Common evaluation time of the fields for single time frame.
    std::vector<double> observe_values_time_;

    // Name of the observation stream. Base for the output filename.
    std::string observe_name_;
    
    /// Output file stream.
    std::ofstream observe_file_;

    /// String representation of the time unit.
    std::string time_unit_str_;
    /// Time unit in seconds.
    double time_unit_seconds_;
    /// Precision of float output
    unsigned int precision_;
    
    // Warn for no observe fields only once.
    bool no_fields_warning=false;

	/// Index set assigning to local point index its global index.
    std::vector<LongIdx> point_4_loc_;

	/// Parallel distribution of observe points.
	Distribution *point_ds_;

	/// Index of actual (last) time in \p observe_values_time_ vector
	unsigned int observe_time_idx_;

	friend class ObservePointAccessor;
};


/**
 * @brief Point accessor allow iterate over local Observe points.
 *
 * Iterator is defined by:
 *  - Observe object
 *  - local index of observe point (iterated value)
 */
class ObservePointAccessor {
public:
    /// Default invalid accessor.
	ObservePointAccessor()
    : observe_(NULL)
    {}

    /**
     * Observe point accessor.
     */
	ObservePointAccessor(const Observe *observe, unsigned int loc_idx)
    : observe_(observe), loc_point_idx_(loc_idx)
    {}

    /// Return local index to point.
    inline unsigned int local_idx() const {
        return loc_point_idx_;
    }

    /// Return global index to point.
    inline unsigned int global_idx() const {
        return observe_->point_4_loc_[loc_point_idx_];
    }

    /// Return ElementAccessor to element of loc_ele_idx_.
    inline const ObservePoint observe_point() const {
    	return observe_->points_[ this->global_idx() ];
    }

    /// Return local index in data cache (combination of local point index and index of stored time)
    inline unsigned int loc_point_time_index() const {
        return (observe_->point_4_loc_.size() * observe_->observe_time_idx_) + loc_point_idx_;
    }

    /// Check validity of accessor (see default constructor)
    inline bool is_valid() const {
        return observe_ != NULL;
    }

    /// Iterates to next local point.
    inline void inc() {
    	loc_point_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const ObservePointAccessor& other) {
    	return (loc_point_idx_ == other.loc_point_idx_);
    }

protected:
    /// Pointer to the Observe owning the point.
    const Observe * observe_;
    /// Index into Observe::point_4_loc_ array.
    unsigned int loc_point_idx_;

};




#endif /* SRC_IO_OBSERVE_HH_ */
