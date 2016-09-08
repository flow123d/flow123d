/*
 * output_time_set.cc
 *
 *  Created on: Jul 11, 2016
 *      Author: jb
 */

#include "tools/time_marks.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "io/output_time_set.hh"
#include "tools/time_governor.hh"

namespace IT = Input::Type;


const IT::Array OutputTimeSet::get_input_type()
{
    static const IT::Record &time_grid =
        IT::Record("TimeGrid", "Equally spaced grid of time points.")
            .allow_auto_conversion("begin")
            .declare_key("begin", IT::Double(0.0), IT::Default::read_time("The initial time of the associated equation."),
                    "The start time of the grid.")
            .declare_key("step", IT::Double(0.0), IT::Default::optional(),
                    "The step of the grid. If not specified, the grid consists only of the start time.")
            .declare_key("end", IT::Double(0.0), IT::Default::read_time("The end time of the simulation."),
                    "The time greater or equal to the last time in the grid.")
            .close();
    return IT::Array(time_grid);
}


void OutputTimeSet::read_from_input(Input::Array in_array, const TimeGovernor &tg)
{
    read_from_input(in_array, tg, tg.equation_fixed_mark_type());
}

void OutputTimeSet::read_from_input(Input::Array in_array, const TimeGovernor &tg, TimeMark::Type mark_type)
{
    double initial_time = tg.init_time();
    double simulation_end_time = tg.end_time();

    for(auto it =in_array.begin<Input::Record>(); it != in_array.end(); ++it) {
        double t_begin = it->val<double>("begin", initial_time);
        double t_end = it->val<double>("end", simulation_end_time );
        double t_step;
        if (! it->opt_val("step", t_step) ) {
            t_end = t_begin;
            t_step = 1.0;
        }
        if ( t_begin > t_end) {
            WarningOut().fmt("Ignoring output time grid. Time begin {}  > time end {}. {}",
                    t_begin, t_end, it->address_string());
            continue;
        }
        if ( t_step < 2*numeric_limits<double>::epsilon()) {
            WarningOut().fmt("Ignoring output time grid. Time step {} < two times machine epsilon {}. {}",
                    t_step, 2*numeric_limits<double>::epsilon(),it->address_string());
            continue;
        }

        this->add(t_begin, t_step, t_end, mark_type);
    }
}


bool OutputTimeSet::contains(TimeMark mark) const {
    return times_.find( mark.time() ) != times_.end();
}


void OutputTimeSet::add(double time, TimeMark::Type mark_type)
{
    TimeMark::Type output_mark_type = mark_type | TimeGovernor::marks().type_output();
    auto mark = TimeMark(time, output_mark_type);
    double mark_time = TimeGovernor::marks().add(mark).time();
    times_.insert( mark_time );
}



void OutputTimeSet::add(double begin, double step, double end, TimeMark::Type mark_type)
{
    //DebugOut().fmt("set add: {} {} {} {}\n", begin, step, end, mark_type);
    ASSERT_GE( step, 2*numeric_limits<double>::epsilon());
    ASSERT_LE( begin, end);
    ASSERT_LT( end, TimeGovernor::inf_time );
    TimeMark::Type output_mark_type = mark_type | TimeGovernor::marks().type_output();

    unsigned int n_steps=((end - begin) / step + TimeGovernor::time_step_precision);
    for(unsigned int i = 0; i <= n_steps; i++) {
        auto mark = TimeMark(begin + i * step, output_mark_type);
        double time = TimeGovernor::marks().add(mark).time();
        times_.insert( time );
        //DebugOut().fmt("add time: {} size: {}\n", time, times_.size());
    }
}
