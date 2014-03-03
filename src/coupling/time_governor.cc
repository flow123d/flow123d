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

//initialize constant pointer to TimeMarks object
TimeMarks * const TimeGovernor::time_marks = new TimeMarks();

// fraction of subsequent time steps can not be less then the comparison_precision
const double TimeGovernor::comparison_precision = 0.0001;
const double TimeGovernor::time_step_lower_bound = numeric_limits<double>::epsilon();
const double TimeGovernor::inf_time =  numeric_limits<double>::infinity();
const double TimeGovernor::round_n_steps_precision = 1e-14;



using namespace Input::Type;

Record TimeGovernor::input_type = Record("TimeGovernor",
            "Setting of the simulation time. (can be specific to one eqaution)")
    .declare_key("start_time", Double(), Default("0.0"),
                "Start time of the simulation.")
    .declare_key("end_time", Double(), Default::obligatory(),
                "End time of the simulation.")
    .declare_key("init_dt", Double(0.0), Default::optional(),
                                    "Initial guess for the time step. The time step is fixed if "
                                    "hard time step limits are not set.")
    .declare_key("min_dt", Double(0.0), Default::read_time("Machine precision or 'init_dt' if specified"),
                                    "Hard lower limit for the time step.")
    .declare_key("max_dt", Double(0.0), Default::read_time("Whole time of the simulation or 'init_dt' if specified"),
                                    "Hard upper limit for the time step.");




/*
 * TODO:
 * TimeGovernor should be constructed from JSON object.
 */
//TimeGovernor::TimeGovernor(const Input::Record &input, TimeMarks &marks, const TimeMark::Type
TimeGovernor::TimeGovernor(const Input::Record &input, const TimeMark::Type fixed_time_mask)
: time_level(0),
  time_step(time_step_lower_bound),
  last_time_step(time_step_lower_bound),
  fixed_dt(0.0),
  dt_fixed_now(false),
  dt_changed(true),
  upper_constraint_(inf_time),
  lower_constraint_(time_step_lower_bound),
  max_time_step(inf_time),
  min_time_step(time_step_lower_bound),
  //time_marks(&marks),
  fixed_time_mark_mask(fixed_time_mask | time_marks->type_fixed_time()),
  steady(false)
{
    time = input.val<double>("start_time");
    if (time < 0.0) xprintf(UsrErr, "Start time has to be equal or greater than ZERO.\n");
    
    end_time_ = input.val<double>("end_time");
    if (end_time_ < time) xprintf(UsrErr, "End time must be greater than start time.\n");
    
    end_of_fixed_dt_interval = time;
    //time_step_constraint = end_time_ - time;

    if (end_time_ != inf_time)  max_time_step = end_time_ - time;

    time_step = max_time_step;

    time_marks->add( TimeMark(time, fixed_time_mark_mask) );
    time_marks->add( TimeMark(end_time_, fixed_time_mark_mask) );

    Input::Iterator<double> it = input.find<double>("init_dt");
    if (it) {
        if (*it <= 0.0) xprintf(UsrErr, "Initial time step has to be greater than ZERO.\n");
        time_step = *it;
        max_time_step = time_step;
        min_time_step = time_step;
        upper_constraint_ = time_step;
        lower_constraint_ = time_step;
    }
    
    max_time_step = input.val<double>("max_dt", max_time_step);
    min_time_step = input.val<double>("min_dt", min_time_step);

    if (min_time_step <= 0.0) xprintf(UsrErr, "Minimal time step has to be greater than ZERO.\n");
    if (max_time_step < min_time_step) xprintf(UsrErr, "Maximal time step has to be greater or equal to the minimal.\n");
    
    last_time_step=0.0;
}


TimeGovernor::TimeGovernor(double init_time, double dt)
: time_level(0),
  time_step(dt),
  last_time_step(time_step_lower_bound),
  fixed_dt(dt),
  dt_fixed_now(true),
  dt_changed(true),
  upper_constraint_(inf_time),
  lower_constraint_(time_step_lower_bound),
  max_time_step(dt),
  min_time_step(dt),
  //time_marks(&marks),
  fixed_time_mark_mask(time_marks->type_fixed_time()),
  steady(false)
{
    time = init_time;
    if (time < 0.0) xprintf(UsrErr, "Start time has to be equal or greater than ZERO.\n");

    end_time_ = inf_time;
    end_of_fixed_dt_interval = inf_time;

    if (min_time_step <= 0.0) xprintf(UsrErr, "Minimal time step has to be greater than ZERO.\n");
    if (max_time_step < min_time_step) xprintf(UsrErr, "Maximal time step has to be greater or equal to the minimal.\n");

    last_time_step=0.0;
}


// steady time governor constructor
TimeGovernor::TimeGovernor()
: time_level(0),
  time(0.0),
  end_of_fixed_dt_interval(time),
  end_time_(inf_time),
  time_step(inf_time),
  last_time_step(0.0),
  fixed_dt(0.0),
  dt_changed(true),
  upper_constraint_(inf_time),
  lower_constraint_(time_step_lower_bound),
  max_time_step(inf_time),
  min_time_step(time_step_lower_bound),
  fixed_time_mark_mask(0x0),
  steady(true)
{
    time_marks->add( TimeMark(time, fixed_time_mark_mask) );
}


TimeGovernor::TimeGovernor(double init_time)
: time_level(0),
  time(init_time),
  end_of_fixed_dt_interval(time),
  end_time_(inf_time),
  time_step(inf_time),
  last_time_step(0.0),
  fixed_dt(0.0),
  dt_changed(true),
  upper_constraint_(inf_time),
  lower_constraint_(time_step_lower_bound),
  max_time_step(inf_time),
  min_time_step(time_step_lower_bound),
  fixed_time_mark_mask(0x0),
  steady(true)
{
    time_marks->add( TimeMark(init_time, fixed_time_mark_mask) );
}




void TimeGovernor::set_permanent_constraint( double min_dt, double max_dt)
{
    if (min_dt < 0.0) xprintf(UsrErr, "Minimal time step has to be greater than ZERO.\n");
    if (max_dt < min_dt) xprintf(UsrErr, "Maximal time step has to be greater or equal to the minimal.\n");

    min_time_step = max(min_dt, time_step_lower_bound);
    max_time_step = min(max_dt, end_time_-time);
    upper_constraint_ = max_time_step;
    lower_constraint_ = min_time_step;
}


// int set_constain - dle nastaveni constraint
// interval - constraint - jako v cmp u Qsortu
// -1 vetsi nez interval (min_dt,max_dt)
// +1 mensi
// 0 OK

int TimeGovernor::set_upper_constraint (double upper)
{
    if (upper_constraint_ < upper) 
    {
        //do not change upper_constraint_
        return -1;
    }
    
    if (lower_constraint_ <= upper) 
    {
        //change upper_constraint_ to upper
        upper_constraint_ = upper;
        return 0;
    }
    
    if (lower_constraint_ > upper) 
    {
        //do not change upper_constraint_
        return 1;
    }

    return 0;
}

int TimeGovernor::set_lower_constraint (double lower)
{   
    if (upper_constraint_ < lower) 
    {
        //do not change lower_constraint_
        return -1;
    }
    
    if (min_time_step <= lower) 
    {
        //change lower_constraint_ to lower
        lower_constraint_ = lower;
        return 0;
    }
    
    if (min_time_step > lower) 
    {
        //do not change lower_constraint_
        return 1;
    }

    return 0;
}


double TimeGovernor::estimate_dt() const {
    if (time == inf_time || is_end()) return 0.0;
    if (time_marks == NULL) return inf_time;
    
    //two states of steady time governor returns different estimate
    if (steady && time_level==0) return inf_time;	//at the beginning
    if (steady && time_level==1) return 0.0;		//at the end
    
    if (this->lt(end_of_fixed_dt_interval))    return fixed_dt;

    // jump to the first future fix time
    TimeMarks::iterator fix_time_it = time_marks->next(*this, fixed_time_mark_mask);
    // compute step to next fix time and apply constraints
    double full_step = fix_time_it->time() - time;
    
    double step_estimate = min(full_step, upper_constraint_);
    step_estimate = max(step_estimate, lower_constraint_); //these two must be in this order

    // round the time step to have integer number of steps till next fix time
    // this always select shorter time step
    int n_steps = ceil( full_step / step_estimate );
    //xprintf(Msg, "float_n_step: %.16f ", full_step / step_estimate);
    
    int n_floor_steps = floor(full_step / step_estimate);
    
    //checking rounding precision: e.g. 1.0000000000000009 -> 1 NOT 2
    if ( abs(full_step / step_estimate - n_floor_steps) < round_n_steps_precision)
	n_steps = n_floor_steps;
    
    step_estimate = full_step / n_steps;
    //xprintf(Msg, "n_step: %d step_estimate: %f\n", n_steps, step_estimate);
    
    // if the step estimate gets by rounding lower than lower constraint program will not crash
    // will just output a user warning.
    if (step_estimate < lower_constraint_)
        xprintf(Warn, "Time step estimate is below the lower constraint of time step. The difference is: %.16f.\n", 
                lower_constraint_ - step_estimate);
    
    /*    OBSOLETE
    // check permanent bounds with a bit of tolerance
    if (step_estimate < min_time_step*0.99) {
        // try longer step
        double longer_step = full_step / (n_steps - 1);
        if (longer_step <= max_time_step*1.01) {
            step_estimate = longer_step;
        }
    }
    */ 
    return step_estimate;
}



void TimeGovernor::next_time()
{
    if (time == inf_time || is_end()) return;
    
    //in case the time governor is steady the time is set to end time which is infinity
    if (steady) {
        last_time_step = end_time_;
        time_level = 1;
	time_step = 0.0;
        time = end_time_;
        dt_changed = false;
        return;
    }
    
    if (this->lt(end_of_fixed_dt_interval)) {
        // this is done for fixed step
        // make tiny correction of time step in order to avoid big rounding errors
        // tiny correction means that dt_changed 'is NOT changed'
        fixed_dt = (end_of_fixed_dt_interval-time) / round( (end_of_fixed_dt_interval-time) / fixed_dt );
        last_time_step = time_step;
        time_step = fixed_dt;
        
        //checking whether fixed time step has been changed (by fix_dt_until_mark() method) since last time
        if (dt_fixed_now)
        {
            dt_fixed_now = false;       
            
            //is true unless new fixed_dt is not equal previous time_step
            dt_changed = (last_time_step != time_step); 
        }
        else
            dt_changed = false;
    }
    else
    {
        // this is done if end_of_fixed_dt_interval is not set (means it is equal to -infinity)
        last_time_step = time_step;
        time_step = estimate_dt();
        dt_changed= (last_time_step != time_step);
    }
    
    last_time_ = time;
    time+=time_step;
    time_level++;
    // refreshing the upper_constraint_
    upper_constraint_ = min(end_time_ - time, max_time_step);
    lower_constraint_ = min_time_step;
}

