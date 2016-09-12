/*
 * equation_output.hh
 *
 *  Created on: Jul 8, 2016
 *      Author: jb
 */

#ifndef SRC_IO_EQUATION_OUTPUT_HH_
#define SRC_IO_EQUATION_OUTPUT_HH_

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


private:
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


    /// output stream (may be shared by more equation)
    std::shared_ptr<OutputTime> stream_;
    /// The time mark type of the equation.
    TimeMark::Type equation_type_;
    /// The fixed time mark type of the equation.
    TimeMark::Type equation_fixed_type_;
    /// The time set used for the fields without explicit time set.
    OutputTimeSet common_output_times_;

    /// Time sets of individual fields.
    std::unordered_map<string, OutputTimeSet> field_output_times_;

    /// Set of observed fields. The observe points are given within the observe stream.
    std::unordered_set<string> observe_fields_;
};


#endif /* SRC_IO_EQUATION_OUTPUT_HH_ */
