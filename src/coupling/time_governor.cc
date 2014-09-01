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
TimeMarks * const TimeGovernor::time_marks_ = new TimeMarks();

// fraction of subsequent time steps can not be less then the comparison_precision
const double TimeGovernor::comparison_precision = 0.0001;
const double TimeGovernor::time_step_lower_bound = numeric_limits<double>::epsilon();
const double TimeGovernor::inf_time =  numeric_limits<double>::infinity();
const double TimeGovernor::round_n_steps_precision = 1e-14;



using namespace Input::Type;





Record TimeGovernor::input_type = Record("TimeGovernor",
            "Setting of the simulation time. (can be specific to one equation)")
	.allow_auto_conversion("init_dt")
    .declare_key("start_time", Double(), Default("0.0"),
                "Start time of the simulation.")
    .declare_key("end_time", Double(), Default::optional(),
                "End time of the simulation.")
    .declare_key("init_dt", Double(0.0), Default("0.0"),
            "Initial guess for the time step.\n"
    		"The time step is fixed if the hard time step limits are not set.\n"
            "If set to 0.0, the time step is determined in fully autonomous way if the equation supports it.")
    .declare_key("min_dt", Double(0.0), Default::read_time("Machine precision or 'init_dt' if specified"),
                                    "Hard lower limit for the time step.")
    .declare_key("max_dt", Double(0.0), Default::read_time("Whole time of the simulation or 'init_dt' if specified"),
                                    "Hard upper limit for the time step.");




/*
 * TODO:
 * TimeGovernor should be constructed from JSON object.
 */
//TimeGovernor::TimeGovernor(const Input::Record &input, TimeMarks &marks, const TimeMark::Type
TimeGovernor::TimeGovernor(const Input::Record &input, TimeMark::Type eq_mark_type)
{
	// use INF end time as default
	end_time_ = inf_time;
    if (input.opt_val<double>("end_time", end_time_))

    // use new mark type as default
    if (eq_mark_type == TimeMark::none_type) eq_mark_type = marks().new_mark_type();

    try {
    	init_common(input.val<double>("init_dt"),
    				input.val<double>("start_time"),
    				end_time_,
    				eq_mark_type);

    	// possibly overwrite limits
    	set_permanent_constraint(
    		input.val<double>("min_dt", min_time_step_),
    		input.val<double>("max_dt", max_time_step_)
    		);
    } catch(ExcTimeGovernorMessage &exc) {
    	exc << input.ei_address();
    	throw;
    }

}

TimeGovernor::TimeGovernor(double init_time, double dt)
{
	init_common(dt, init_time, inf_time, TimeMark::none_type);
}


// steady time governor constructor
TimeGovernor::TimeGovernor(double init_time, TimeMark::Type eq_mark_type)
{
    // use new mark type as default
    if (eq_mark_type == TimeMark::none_type) eq_mark_type = marks().new_mark_type();

	init_common(0.0, init_time, inf_time, eq_mark_type);
	steady_ = true;
}



// common part of constructors
void TimeGovernor::init_common(double dt, double init_time, double end_time, TimeMark::Type type)
{
	time_level_=0;

    if (init_time < 0.0) {
		THROW(ExcTimeGovernorMessage()
				<< EI_Message("Start time has to be greater or equal to 0.0\n")

				);
    }

	init_time_  = time_ = init_time;
	last_time_ = -inf_time;
	last_time_step_ = inf_time;


	if (end_time < init_time) {
		THROW(ExcTimeGovernorMessage()	<< EI_Message("End time must be greater than start time.\n") );
    }

    end_time_ = end_time;

    if (dt == 0.0) {
    	// variable time step
    	fixed_time_step_=0.0;
    	is_time_step_fixed_=false;
    	time_step_changed_=true;
    	end_of_fixed_dt_interval_ = time_;

    	min_time_step_=lower_constraint_=time_step_lower_bound;
    	if (end_time_ == inf_time) {
        	max_time_step_=upper_constraint_=inf_time;
    	} else {
    		max_time_step_=upper_constraint_= end_time - time_;
    	}
    	// choose maximum possible time step
    	time_step_=max_time_step_;
    } else {
    	// fixed time step
    	if (dt < time_step_lower_bound)
    		xprintf(UsrErr, "Fixed time step small then machine precision. \n");

    	fixed_time_step_=dt;
    	is_time_step_fixed_=true;
    	time_step_changed_=true;
    	end_of_fixed_dt_interval_ = inf_time;

    	upper_constraint_=max_time_step_=dt;
    	lower_constraint_=min_time_step_=dt;
    	time_step_=dt;

    }

	eq_mark_type_=type;
	steady_=false;
    time_marks_->add( TimeMark(time_, equation_fixed_mark_type()) );
    if (end_time_ != inf_time)
    	time_marks_->add( TimeMark(end_time_, equation_fixed_mark_type()) );
}






void TimeGovernor::set_permanent_constraint( double min_dt, double max_dt)
{
    if (min_dt < time_step_lower_bound) {
		THROW(ExcTimeGovernorMessage()	<< EI_Message("'min_dt' smaller then machine precision.\n") );
    }
    if (max_dt < min_dt) {
		THROW(ExcTimeGovernorMessage()	<< EI_Message("'max_dt' smaller then 'min_dt'.\n") );
    }

    lower_constraint_ = min_time_step_ = max(min_dt, time_step_lower_bound);
    upper_constraint_ = max_time_step_ = min(max_dt, end_time_-time_);
}


// int set_constrain - dle nastaveni constraint
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
    
    if (min_time_step_ <= lower) 
    {
        //change lower_constraint_ to lower
        lower_constraint_ = lower;
        return 0;
    }
    
    if (min_time_step_ > lower) 
    {
        //do not change lower_constraint_
        return 1;
    }

    return 0;
}



void TimeGovernor::add_time_marks_grid(double step, TimeMark::Type mark_type) const
{
	if (end_time() == inf_time) {
		THROW(ExcTimeGovernorMessage()
				<< EI_Message("Missing end time for make output grid required by key 'time_step' of the output stream.\n")
			);
	}

	marks().add_time_marks(init_time_, step, end_time(), mark_type | eq_mark_type_);
	// always add start time and end time
	marks().add(TimeMark(init_time_, mark_type | eq_mark_type_));
	marks().add(TimeMark(end_time(), mark_type | eq_mark_type_));
}



double TimeGovernor::estimate_dt() const {
    if (time_ == inf_time || is_end()) return 0.0;
    if (time_marks_ == NULL) return inf_time;
    
    //two states of steady time governor returns different estimate
    //if (steady && time_level==0) return inf_time;	//at the beginning
    //if (steady && time_level==1) return 0.0;		//at the end
    
    if (this->lt(end_of_fixed_dt_interval_))    return fixed_time_step_;

    // jump to the first future fix time
    TimeMarks::iterator fix_time_it = time_marks_->next(*this, equation_fixed_mark_type());
    // compute step to next fix time and apply constraints
    double full_step = fix_time_it->time() - time_;
    
    double step_estimate = min(full_step, upper_constraint_);
    step_estimate = max(step_estimate, lower_constraint_); //these two must be in this order

    if (step_estimate == inf_time) return step_estimate;

    // round the time step to have integer number of steps till next fix time
    // this always select shorter time step
    int n_steps = ceil( full_step / step_estimate );
    //xprintf(Msg, "float_n_step: %.16f ", full_step / step_estimate);
    
    
    //checking rounding precision: e.g. 1.0000000000000009 -> 1 NOT 2
    int n_floor_steps = floor(full_step / step_estimate);
    if ( abs(full_step / step_estimate - n_floor_steps) < round_n_steps_precision)   	n_steps = n_floor_steps;
    
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
    ASSERT_LE(0.0, time_);
    if (time_ == inf_time || is_end()) return;
    
    //in case the time governor is steady the time is set to end time which is infinity
   /* if (steady) {
        last_time_ = time;
        last_time_step = end_time_;
        time_level = 1;
        time_step = 0.0;
        time = end_time_;
        dt_changed = false;
        return;
    }*/
    
    if (this->lt(end_of_fixed_dt_interval_)) {
        // this is done for fixed step
        // make tiny correction of time step in order to avoid big rounding errors
        // tiny correction means that dt_changed 'is NOT changed'
    	if (end_of_fixed_dt_interval_ < inf_time) {
    		fixed_time_step_ = (end_of_fixed_dt_interval_-time_) / round( (end_of_fixed_dt_interval_-time_) / fixed_time_step_ );
    	}

    	last_time_step_ = time_step_;
        time_step_ = fixed_time_step_;
        
        //checking whether fixed time step has been changed (by fix_dt_until_mark() method) since last time
        if (is_time_step_fixed_)
        {
            is_time_step_fixed_ = false;       
            
            //is true unless new fixed_dt is not equal previous time_step
            time_step_changed_ = (last_time_step_ != time_step_); 
        }
        else
            time_step_changed_ = false;
    }
    else
    {
        // this is done if end_of_fixed_dt_interval is not set (means it is equal to -infinity)
        last_time_step_ = time_step_;
        time_step_ = estimate_dt();
        time_step_changed_= (last_time_step_ != time_step_);
    }

    last_time_ = time_;
    time_+=time_step_;
    time_level_++;
    // refreshing the upper_constraint_
    upper_constraint_ = min(end_time_ - time_, max_time_step_);
    lower_constraint_ = min_time_step_;

    //if (time_level_ == 7) time_ = end_time_;
}



void TimeGovernor::view(const char *name) const
{
    //xprintf(Msg, "\nTG[%s]: level: %d end_time: %f time: %f step: %f upper: %f lower: %f end_fixed_time: %f type: %x\n",
    //        name, time_level_, end_time_, time_, time_step_, upper_constraint_, lower_constraint_, end_of_fixed_dt_interval_, eq_mark_type_);

    xprintf(Msg, "\nTG[%s]:%06d    t:%10.4f    dt:%10.6f    dt_int<%10.6f,%10.6f>",
            name, time_level_, time_, time_step_, lower_constraint_, upper_constraint_ );
#ifdef DEBUG_MESSAGES
    xprintf(Msg, "    end_time: %f end_fixed_time: %f type: 0x%x\n" , end_time_,  end_of_fixed_dt_interval_, eq_mark_type_);
#else
    xprintf(Msg,"\n");
#endif
}
