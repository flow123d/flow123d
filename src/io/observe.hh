/*
 * observe.hh
 *
 *  Created on: Jun 28, 2016
 *      Author: jb
 */

#ifndef SRC_IO_OBSERVE_HH_
#define SRC_IO_OBSERVE_HH_

#include <string>

#include "input/input_type.hh"
#include "fields/field.hh"
#include "io/output_data.hh"
#include <armadillo>



/**
 * Internal class for the class Observe.
 * Members: input_pos_, snap_dim_, snap_region_name_ are set in constructor. Should be checked before passed in.
 * Members: element_idx_, global_coords_, local_coords_ are derived, set in Observe::find_observe_points.
 */
class ObservePoint {
public:
    DECLARE_INPUT_EXCEPTION(ExcNoInitialPoint,
            << "Failed to find the element containing the initial observe point.\n");
    static const Input::Type::Record & get_input_type();

protected:
    /// Default constructor just for testing.
    ObservePoint();

    /**
     * Constructor. Read from input.
     */
    ObservePoint(Input::Record in_rec, unsigned int point_idx);

    /**
     * Update the observe element and the projection of the initial point on it.
     */
    void update_projection(unsigned int i_elm, arma::vec local_coords, arma::vec3 global_coords);

    /**
     * Returns true if we have already found any observe element.
     */
    bool have_observe_element();

    /**
     * Snap local coords to the subelement. Called by the snap method.
     */
    template <int ele_dim>
    void snap_to_subelement();

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


    /// Index in the input array.
    Input::Record in_rec_;

    /// Observation point name.
    std::string name_;

    /// Input coordinates of the initial position of the observation point.
    arma::vec3 input_point_;

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
     * Maximal number of observe element search levels.
     */
    unsigned int max_levels_;

    /// Final element of the observe point. The index in the mesh.
    unsigned int element_idx_;

    /// Global coordinates of the observation point.
    arma::vec3 global_coords_;

    /// Local (barycentric) coordinates on the element.
    arma::vec local_coords_;

    /// Distance of found projection from the initial point.
    /// If we find more candidates we pass in the closest one.
    double distance_;



    /// Only Observe should use this class directly.
    friend class Observe;

};




class Observe {
public:

    /**
     * Construct the observation object.
     *
     * observe_name - base name of the output file.
     * mesh -
     */
    Observe(string observe_name, Mesh &mesh, Input::Record in_rec);

    /// Destructor, must close the file.
    ~Observe();


    /**
     * Evaluates and store values of the given field in the observe points.
     */
    template<int spacedim, class Value>
    void compute_field_values(Field<spacedim, Value> &field);

    /**
     * Provides a vector of element indices on which the observation values are evaluated.
     * This can be used to evaluate derived fields only on these elements in the times not selected to
     * full output.
     */
    inline const std::set<unsigned int> &observed_elements() const
            { return observed_element_indices_;}

    /**
     * Write field values to the output file. Using the YAML format.
     */
    void output_time_frame();



protected:
    Mesh *mesh_;
    Input::Record in_rec_;


    /// Full information about observe points.
    std::vector<ObservePoint> points_;
    /// Elements of the o_points.
    std::set<unsigned int> observed_element_indices_;

    typedef std::shared_ptr<OutputDataBase> OutputDataPtr;
    typedef std::map< string,  OutputDataPtr > OutputDataFieldMap;

    /// Stored field values.
    OutputDataFieldMap observe_field_values_;

    /// Fields to observe.
    std::set<string> field_names_;

    /// Output file stream.
    std::ofstream observe_file_;



};



#endif /* SRC_IO_OBSERVE_HH_ */
