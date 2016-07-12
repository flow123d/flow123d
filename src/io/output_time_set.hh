/*
 * output_time_set.hh
 *
 *  Created on: Jul 11, 2016
 *      Author: jb
 */

#ifndef SRC_IO_OUTPUT_TIME_SET_HH_
#define SRC_IO_OUTPUT_TIME_SET_HH_


//#include  "boost/unordered_set.hpp"
#include "tools/time_marks.hh"
#include <set>

class OutputTimeSet {
public:

    //OutputTimeSet &operator=(const OutputTimeSet &);

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
    bool contains(TimeMark mark) const;

    void add(double begin, double step, double end, TimeMark::Type mark_type);



private:
    //boost::unordered::unordered_set<TimeMark, TimeMarkHash> time_marks_;
    std::set<double> times_;
};





#endif /* SRC_IO_OUTPUT_TIME_SET_HH_ */
