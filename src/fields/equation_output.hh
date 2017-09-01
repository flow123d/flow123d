/*
 * equation_output.hh
 *
 *  Created on: Jul 8, 2016
 *      Author: jb
 */

#ifndef SRC_FIELDS_EQUATION_OUTPUT_HH_
#define SRC_FIELDS_EQUATION_OUTPUT_HH_

#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "tools/time_marks.hh"
#include "fields/field_set.hh"
#include "input/input_type_forward.hh"
#include "input/accessors_forward.hh"
#include "io/output_time_set.hh"
class TimeStep;


/**
 * A class  responsible for check for output times of individual fields
 * and store their values into the connected output stream.
 */
class EquationOutput : public FieldSet {
public:

    DECLARE_EXCEPTION(ExcFieldNotScalar, << "Field '" << FieldCommon::EI_Field::qval
                                         << "' is not scalar in spacedim 3.");

    /**
     * Make Input::Type for the output record. Particular selection of output fields is created
     * from the contents of *this FeildSet using provided equation name and additional description.
     */
    const Input::Type::Instance &make_output_type(const string &equation_name, const string &aditional_description = "");

    /**
     * Setup the object. Set output stream for field and observe output, input record for configuration of the object and
     * TimeGovernor. The time governor is used to get the equation time mark type, the initial and the end time of the equation.
     */
    void initialize(std::shared_ptr<OutputTime> stream, Input::Record in_rec, const TimeGovernor & tg);

    /**
     * Returns true if @param field is marked for output in the given time @param step.
     */
    bool is_field_output_time(const FieldCommon &field, TimeStep step) const;

    /**
     * Performs output of the fields marked for output in the time @param step.
     */
    void output(TimeStep step);

    /// Selects the error control field out of output field set according to input record.
    void select_error_control_field(std::string error_control_field_name);


private:
    /// Store pair of OutputTimeSet and DiscreteSpace data of Field
    struct OutputTimeData {
    	/// Empty constructor
    	OutputTimeData() {}
    	/// Constructor
    	OutputTimeData(OutputTimeSet times, OutputTime::DiscreteSpace disc)
    	: output_times(times), discrete(disc) {}

    	OutputTimeSet output_times;          ///< Time set of Field
    	OutputTime::DiscreteSpace discrete;  ///< Discrete space (can be set by user)
    };

    /**
     * Input type of the configuration record.
     */
    static Input::Type::Record &get_input_type();

    /**
     * Read from the input, set output times and time marks. Must be called after set_stream.
     * TODO: add output_stream times. Optional or always?
     */
    void read_from_input(Input::Record in_rec, const TimeGovernor & tg);

    /**
     * Add a time grid to the common_output_times.
     */
    void add_output_times(double begin, double step, double end);


    /**
     * Create the output mesh of \p stream_ OutputTime object. The field set passed in is used
     * to select the field used for adaptivity of the output mesh.
     */
    void make_output_mesh();


    /**
     * Get DiscreteSpace type of given field.
     */
    OutputTime::DiscreteSpace get_field_discrete_space(const FieldCommon &field) const;


    /// output stream (may be shared by more equation)
    std::shared_ptr<OutputTime> stream_;
    /// The time mark type of the equation.
    TimeMark::Type equation_type_;
    /// The fixed time mark type of the equation.
    TimeMark::Type equation_fixed_type_;
    /// The time set used for the fields without explicit time set.
    OutputTimeSet common_output_times_;

    /// Time sets of individual fields.
    std::unordered_map<string, OutputTimeData> field_output_times_;

    /// Set of observed fields. The observe points are given within the observe stream.
    std::unordered_set<string> observe_fields_;

    /// Refinement error control field.
    Field<3, FieldValue<3>::Scalar> *error_control_field_;
};


#endif /* SRC_FIELDS_EQUATION_OUTPUT_HH_ */
