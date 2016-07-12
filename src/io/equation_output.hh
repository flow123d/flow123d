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


    void initialize(std::shared_ptr<OutputTime> stream, Input::Record in_rec, const TimeGovernor & tg);
    /**
     * Collective interface to @p FieldCommonBase::output_type().
     * @param rt   Discrete function space (element, node or corner data).
     */
    void set_output_type(OutputTime::DiscreteSpace rt);
    /*{
        for(FieldCommon *field : field_list) field->output_type(rt);
    }*/


    bool is_field_output_time(const FieldCommon &field, TimeStep step) const;

    void output(TimeStep step);


private:
    static Input::Type::Record &get_input_type();

    /**
     * Read from the input, set output times and time marks. Must be called after set_stream.
     * TODO: add output_stream times. Optional or always?
     */
    void read_from_input(Input::Record in_rec, const TimeGovernor & tg);

    //void add_output_time(double begin);
    void add_output_times(double begin, double step, double end);


    std::shared_ptr<OutputTime> stream_;
    TimeMark::Type equation_type_;
    TimeMark::Type equation_fixed_type_;
    OutputTimeSet common_output_times_;
    std::unordered_map<string, OutputTimeSet> field_output_times_;
    std::unordered_set<string> observe_fields_;
};


#endif /* SRC_IO_EQUATION_OUTPUT_HH_ */
