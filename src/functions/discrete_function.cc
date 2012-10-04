/*
 * discrete_function.cc
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */

#include "input/input_type.hh"
#include "input/accessors.hh"
#include "functions/discrete_function.hh"

Input::Type::Record &DiscreteFunction::get_input_type() {
using namespace Input::Type;
    static Record rec("DiscreteFunction",
            "Function that returns natural integers. Currently, the function is either constant or table (i.e. its domain is interval of integral indices).");

    if (! rec.is_finished()) {
        rec.declare_key("constant", Integer(0),
                Default::read_time("If 'constant' is not given 'table' has to be specified."), "Value for constant function.");
        rec.declare_key("table", Array(Integer(0)),
                Default("0"), "Value for table function.");
        rec.allow_auto_conversion("constant");
        rec.finish();
    }
    return rec;
}




DiscreteFunction::DiscreteFunction(const Input::Record &fce_rec) {
    Input::Iterator<unsigned int> it = fce_rec.find<unsigned int>("constant");
    if (it) {
        value_ = *it;
    } else {
        Input::Array vec = fce_rec.val<Input::Array>("table");
        vec.copy_to( values_ );
        value_=0;
    }
}



DiscreteFunction::DiscreteFunction(const std::vector<unsigned int> &values)
: value_(0), values_(values)
{}



unsigned int DiscreteFunction::value(unsigned int idx)
{
    if (values_.size() == 0)
        return value_;
    else {
        ASSERT_LESS( idx, vector_.size());
        return values_[idx];
    }

}


