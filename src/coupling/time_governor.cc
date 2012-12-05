/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id: time_governor.cc 1877 2012-09-27 12:53:36Z jan.brezina $
 * $Revision: 1877 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-09-27 14:53:36 +0200 (Čt, 27 zář 2012) $
 *
 * @file
 * @ingroup application
 * @brief Basic time management class.
 */

#include "system/system.hh"
#include "time_governor.hh"
#include "time_marks.hh"
#include "input/accessors.hh"

#include <limits>

// fraction of subsequent time steps can not be less then the comparison_precision
const double TimeGovernor::comparison_precision = 0.0001;
const double TimeGovernor::time_step_lower_bound = numeric_limits<double>::epsilon();
const double TimeGovernor::inf_time =  numeric_limits<double>::infinity();

/*
 * TODO:
 * TimeGovernor should be constructed from JSON object.
 */
TimeGovernor::TimeGovernor(const Input::Record &input, TimeMarks &marks, const TimeMark::Type fixed_time_mask)
: time_level(0),
  time_step(time_step_lower_bound),
  last_time_step(time_step_lower_bound),
  fixed_dt(0.0),
  dt_changed(true),
  max_time_step(inf_time),
  min_time_step(time_step_lower_bound),
  time_marks(&marks),
  fixed_time_mark_mask(fixed_time_mask | time_marks->type_fixed_time())
{
    time  =  input.val<double>("start_time");
    end_of_fixed_dt_interval=time;
    end_time_  =  input.val<double>("end_time");
    time_step_constrain = end_time_ - time;


    if (end_time_ != inf_time)  max_time_step=end_time_ - time;

    time_step=max_time_step;

    time_marks->add( TimeMark(time, fixed_time_mark_mask) );
    time_marks->add( TimeMark(end_time_, fixed_time_mark_mask) );

    Input::Iterator<double> it = input.find<double>("init_dt");
    if (it) {
        time_step = *it;
        max_time_step = time_step;
        min_time_step = time_step;
    }
    max_time_step = input.val<double>("max_dt", max_time_step);
    min_time_step = input.val<double>("min_dt", min_time_step);

    last_time_step=0.0;
}


// this is get directly form darcy flow staedy,
// TODO: think about generalization
TimeGovernor::TimeGovernor(TimeMarks &marks)
: time_level(0),
  time(-1.0),
  end_of_fixed_dt_interval(time),
  end_time_(inf_time),
  time_step(inf_time),
  last_time_step(0.0),
  fixed_dt(0.0),
  dt_changed(true),
  time_step_constrain(inf_time),
  max_time_step(inf_time),
  min_time_step(time_step_lower_bound),
  time_marks(& marks),
  fixed_time_mark_mask(0x0)

{}



TimeGovernor::TimeGovernor(double init_time)
: time_level(0),
  time(init_time),
  end_of_fixed_dt_interval(time),
  end_time_(inf_time),
  time_step(inf_time),
  last_time_step(0.0),
  fixed_dt(0.0),
  dt_changed(true),
  time_step_constrain(inf_time),
  max_time_step(inf_time),
  min_time_step(time_step_lower_bound),
  time_marks(NULL),
  fixed_time_mark_mask(0x0)

{}


Input::Type::Record &TimeGovernor::get_input_type() {
    using namespace Input::Type;
    static Record rec("TimeGovernor",
            "Setting of the simulation time. (can be specific to one eqaution)");

    if (! rec.is_finished() ) {
        rec.declare_key("start_time", Double(), Default("0.0"),
                "Start time of the simulation.");
        rec.declare_key("end_time", Double(), Default::obligatory(),
                "End time of the simulation.");
        rec.declare_key("init_dt", Double(0.0), Default::optional(),
                                    "Initial guess for the time step. The time step is fixed if "
                                    "hard time step limits are not set.");
        rec.declare_key("min_dt", Double(0.0), Default::read_time("Machine precision or 'init_dt' if specified"),
                                    "Hard lower limit for the time step.");
        rec.declare_key("max_dt", Double(0.0), Default::read_time("Whole time of the simulation or 'init_dt' if specified"),
                                    "Hard upper limit for the time step.");
        rec.finish();
    }
    return rec;
}



void TimeGovernor::set_permanent_constrain( double min_dt, double max_dt)
{
    ASSERT( min_dt >= 0.0,"Minimal time step has to be greater than ZERO\n");
    ASSERT( max_dt >= min_dt,"Maximal time step has to be greater or equal to the minimal.\n");

    min_time_step=max(min_dt, time_step_lower_bound);
    max_time_step=min(max_dt, end_time_-time);
}



void TimeGovernor::set_constrain(double dt_constrain)
{
    time_step_constrain = min(time_step_constrain, dt_constrain);
}



double TimeGovernor::estimate_dt() const {
    if (time == inf_time || is_end()) return 0.0;
    if (time_marks == NULL) return inf_time;
    if (this->lt(end_of_fixed_dt_interval))    return fixed_dt;

    // jump to the first future fix time
    TimeMarks::iterator fix_time_it = time_marks->next(*this, fixed_time_mark_mask);
    // compute step to next fix time and apply constrains
    double full_step = fix_time_it->time() - time;
    double step_estimate = min(full_step, time_step_constrain);
    step_estimate = min(step_estimate, max_time_step);
    step_estimate = max(step_estimate, min_time_step); // possibly overwrites time_step_constrain


    // round the time step to have integer number of steps till next fix time
    // this always select shorter time step
    int n_steps = ceil( full_step / step_estimate );
    step_estimate = full_step / n_steps;

    // check permanent bounds with a bit of tolerance
    if (step_estimate < min_time_step*0.99) {
        // try longer step
        double longer_step = full_step / (n_steps - 1);
        if (longer_step <= max_time_step*1.01) {
            step_estimate = longer_step;
        }
    }

    return step_estimate;
}



void TimeGovernor::next_time()
{
    if (time == inf_time || is_end()) return;
    // TODO: following is because steady solvers but needs better solution
    if (end_time_ == inf_time) {
        time = end_time_;
        return;
    }
    if (this->lt(end_of_fixed_dt_interval)) {
        // make tiny correction of time step in order to avoid big rounding errors
        fixed_dt= end_of_fixed_dt_interval / round( end_of_fixed_dt_interval / fixed_dt );
    }

    //last_time=time;
    last_time_step = time_step;

    time_step = estimate_dt();
    dt_changed= (last_time_step != time_step);
    time_step_constrain = min(end_time_-time, max_time_step); // reset time step constrain

    time+=time_step;
    time_level++;


}

