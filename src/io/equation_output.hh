/*
 * equation_output.hh
 *
 *  Created on: Jul 8, 2016
 *      Author: jb
 */

#ifndef SRC_IO_EQUATION_OUTPUT_HH_
#define SRC_IO_EQUATION_OUTPUT_HH_

#include <memory>
#include  <unordered_map>
#include  "boost/unordered_set.hpp"
#include "tools/time_marks.hh"
#include "fields/field_set.hh"
#include "input/input_type_forward.hh"
#include "input/accessors_forward.hh"

class TimeStep;

class OutputTimeSet {
public:
    /**
     *
     */
    static const Input::Type::Array get_input_type();
    /**
     *
     */
    void read_from_input(Input::Array in_array, TimeMark::Type mark_type);
    /**
     *
     */
    bool contains(TimeMark mark);

private:
    typedef std::pair<double, TimeMark::Type> TimeMarkStorage;

    TimeMark from_time_mark_storage(TimeMarkStorage tm_storage);
    TimeMarkStorage to_time_mark_storage(TimeMark mark);

    TimeMark::Type type_;
    boost::unordered::unordered_set<TimeMarkStorage> time_marks_;
};


/**
 * A class  responsible for check for output times of individual fields
 * and store their values into the connected output stream.
 */
class EquationOutput : public FieldSet {
public:
    EquationOutput(std::shared_ptr<OutputTime> stream, TimeMark::Type equation_mark_type);
    const Input::Type::Record &get_input_type();
    void read_from_input(Input::Record in_rec);
    void output(TimeStep step);

private:
    std::shared_ptr<OutputTime> stream_;
    TimeMark::Type equation_type_;
    OutputTimeSet common_output_times_;
    std::unordered_map<string, OutputTimeSet> field_output_times_;
    boost::unordered::unordered_set<string> observe_fields_;
};


#endif /* SRC_IO_EQUATION_OUTPUT_HH_ */
