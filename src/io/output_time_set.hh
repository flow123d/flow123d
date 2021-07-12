/*
 * output_time_set.hh
 *
 *  Created on: Jul 11, 2016
 *      Author: jb
 */

#ifndef SRC_IO_OUTPUT_TIME_SET_HH_
#define SRC_IO_OUTPUT_TIME_SET_HH_


#include <set>                  // for set
#include "input/type_base.hh"   // for Array
#include "tools/time_marks.hh"  // for TimeMark

class TimeGovernor;
namespace Input { class Array; }


/**
 * Set of times. Simple extension of std::set<double> providing
 * initialization by an array of time grids.
 * TODO: replace by std::set with non-member read and initialize functions.
 *
 * Finally we store just doubles, since it is like a set of the time marks with same type.
 * So we need not to save the type.
 */
class OutputTimeSet {
public:

    /**
     *
     */
    static const Input::Type::Array get_input_type();

    /**
     *
     */
    void read_from_input(Input::Array in_array, const TimeGovernor &tg);
    void read_from_input(Input::Array in_array, const TimeGovernor &tg, TimeMark::Type mark_type);
    /**
     *
     */
    bool contains(TimeMark mark) const;

    void add(double begin, TimeMark::Type mark_type);
    void add(double begin, double step, double end, TimeMark::Type mark_type);



private:
    std::set<double> times_;
};





#endif /* SRC_IO_OUTPUT_TIME_SET_HH_ */
