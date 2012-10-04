/*
 * discrete_function.hh
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */

#ifndef DISCRETE_FUNCTION_HH_
#define DISCRETE_FUNCTION_HH_

#include "input/input_type.hh"
#include "input/accessors.hh"


class DiscreteFunction {
public:
    static Input::Type::Record &get_input_type();

    DiscreteFunction(const Input::Record &fce_rec);
    DiscreteFunction(const std::vector<unsigned int> &values);
    unsigned int value(unsigned int idx);
private:
    unsigned int value_;
    std::vector<unsigned int> values_;
};


#endif /* DISCRETE_FUNCTION_HH_ */
