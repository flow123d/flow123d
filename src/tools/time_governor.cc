/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    time_governor.cc
 * @ingroup application
 * @brief   Basic time management class.
 */

#include  <limits>
#include "system/system.hh"
#include "input/accessors.hh"
#include "time_governor.hh"
#include "time_marks.hh"

//initialize constant pointer to TimeMarks object
TimeMarks TimeGovernor::time_marks_ = TimeMarks();

//const double TimeGovernor::time_step_lower_bound = numeric_limits<double>::epsilon();
const double TimeGovernor::inf_time =  numeric_limits<double>::infinity();
const double TimeGovernor::time_step_precision = 16*numeric_limits<double>::epsilon();


using namespace Input::Type;





const Record & TimeGovernor::get_input_type() {
	return Record("TimeGovernor",
            "Setting of the simulation time. (can be specific to one equation)")
		.allow_auto_conversion("max_dt")
		.declare_key("start_time", Double(), Default("0.0"),
					"Start time of the simulation.")
		.declare_key("end_time", Double(), Default::read_time("Infinite end time."),
					"End time of the simulation.")
		.declare_key("init_dt", Double(0.0), Default("0.0"),
				"Initial guess for the time step.\n"
				"Only useful for equations that use adaptive time stepping."
				"If set to 0.0, the time step is determined in fully autonomous"
				" way if the equation supports it.")
		.declare_key("min_dt", Double(0.0),
				Default::read_time("Machine precision."),
				"Soft lower limit for the time step. Equation using adaptive time stepping can not"
				"suggest smaller time step, but actual time step could be smaller in order to match "
				"prescribed input or output times.")
		.declare_key("max_dt", Double(0.0),
				Default::read_time("Whole time of the simulation if specified, infinity else."),
				"Hard upper limit for the time step. Actual length of the time step is also limited"
				"by input and output times.")
		.close();
}



TimeStep::TimeStep(double init_time) :
index_(0),
length_(1.0),
end_(init_time)
{}



TimeStep::TimeStep() :
index_(0),
length_(TimeGovernor::inf_time),
end_(-TimeGovernor::inf_time)
{}




TimeStep::TimeStep(const TimeStep &other):
index_(other.index_),
length_(other.length_),
end_(other.end_)
{}



TimeStep TimeStep::make_next(double new_length) const
{
    return make_next(new_length, this->end_+new_length);
}



TimeStep TimeStep::make_next(double new_lenght, double end_time) const
{
    TimeStep ts;
    ts.index_=this->index_ +1;
    ts.length_=new_lenght;
    ts.end_=end_time;
    return ts;
}



bool TimeStep::safe_compare(double t1, double t0) const
{
    return t1 >= t0
            - 16*numeric_limits<double>::epsilon()*(1.0+max(abs(t1),abs(t0)));
}



ostream& operator<<(ostream& out, const TimeStep& t_step) {
    out << "time: " << t_step.end() << "step: " << t_step.length() << endl;
}



TimeGovernor::TimeGovernor(const Input::Record &input, TimeMark::Type eq_mark_type)
{
    // use new mark type as default
    if (eq_mark_type == TimeMark::none_type) eq_mark_type = marks().new_mark_type();

    try {
        // set permanent limits
    	init_common(input.val<double>("start_time"),
    				input.val<double>("end_time", inf_time),
    				eq_mark_type);
        set_permanent_constraint(
            input.val<double>("min_dt", min_time_step_),
            input.val<double>("max_dt", max_time_step_)
            );

        double init_dt=input.val<double>("init_dt");
        if (init_dt > 0.0) {
            // set first time step suggested by user
            //time_step_=min(init_dt, time_step_);
            lower_constraint_=init_dt;
            lower_constraint_message_ = "Initial time step set by user.";
            upper_constraint_=init_dt;
            upper_constraint_message_ = "Initial time step set by user.";
        }


    } catch(ExcTimeGovernorMessage &exc) {
    	exc << input.ei_address();
    	throw;
    }

}

TimeGovernor::TimeGovernor(double init_time, double dt)
{
	init_common( init_time, inf_time, TimeMark::none_type);
    // fixed time step
    if (dt < time_step_precision)
        THROW(ExcTimeGovernorMessage() << EI_Message("Fixed time step smaller then machine precision. \n") );

    fixed_time_step_=dt;
    is_time_step_fixed_=true;
    time_step_changed_=true;
    end_of_fixed_dt_interval_ = inf_time;

    lower_constraint_=min_time_step_=dt;
    lower_constraint_message_ = "Initial time step set by user.";
    upper_constraint_=max_time_step_=dt;
    upper_constraint_message_ = "Initial time step set by user.";
    
    //time_step_=dt;
}


// steady time governor constructor
TimeGovernor::TimeGovernor(double init_time, TimeMark::Type eq_mark_type)
{
    // use new mark type as default
    if (eq_mark_type == TimeMark::none_type) eq_mark_type = marks().new_mark_type();

	init_common(init_time, inf_time, eq_mark_type);
	steady_ = true;
}



// common part of constructors
void TimeGovernor::init_common(double init_time, double end_time, TimeMark::Type type)
{


    if (init_time < 0.0) {
		THROW(ExcTimeGovernorMessage()
				<< EI_Message("Start time has to be greater or equal to 0.0\n")

				);
    }

    recent_steps_.set_capacity(size_of_recent_steps_);
    recent_steps_.push_front(TimeStep(init_time));
    init_time_=init_time;

	if (end_time < init_time) {
		THROW(ExcTimeGovernorMessage()	<< EI_Message("End time must be greater than start time.\n") );
    } else {
        end_time_ = end_time;
    }

    //if (dt == 0.0) {
    	// variable time step
    	fixed_time_step_=0.0;
    	is_time_step_fixed_=false;
    	time_step_changed_=true;
    	end_of_fixed_dt_interval_ = init_time_;

    	min_time_step_=lower_constraint_=time_step_precision;
        lower_constraint_message_ = "Permanent minimal constraing, default, time_step_precision.";
    	if (end_time_ == inf_time) {
        	max_time_step_=upper_constraint_=inf_time;
    	} else {
    		max_time_step_=upper_constraint_= end_time - init_time_;
    	}
    	upper_constraint_message_ = "Permanent maximal constraint, default, total simulation time.";
    	// choose maximum possible time step
    	//time_step_=max_time_step_;
    /*} else {
    	// fixed time step
    	if (dt < time_step_lower_bound)
    		THROW(ExcTimeGovernorMessage() << EI_Message("Fixed time step smaller then machine precision. \n") );

    	fixed_time_step_=dt;
    	is_time_step_fixed_=true;
    	time_step_changed_=true;
    	end_of_fixed_dt_interval_ = inf_time;

    	upper_constraint_=max_time_step_=dt;
    	lower_constraint_=min_time_step_=dt;
    	time_step_=dt;

    }*/

	eq_mark_type_=type;
	steady_=false;
    time_marks_.add( TimeMark(init_time_, equation_fixed_mark_type()) );
    if (end_time_ != inf_time)
    	time_marks_.add( TimeMark(end_time_, equation_fixed_mark_type()) );
}






void TimeGovernor::set_permanent_constraint( double min_dt, double max_dt)
{
    if (min_dt < time_step_precision) {
		THROW(ExcTimeGovernorMessage()	<< EI_Message("'min_dt' smaller then machine precision.\n") );
    }
    if (max_dt < min_dt) {
		THROW(ExcTimeGovernorMessage()	<< EI_Message("'max_dt' smaller then 'min_dt'.\n") );
    }

    lower_constraint_ = min_time_step_ = max(min_dt, time_step_precision);
    lower_constraint_message_ = "Permanent minimal constraint, custom.";
    upper_constraint_ = max_time_step_ = min(max_dt, end_time_-t());
    upper_constraint_message_ = "Permanent maximal constraint, custom.";
    
}


// int set_constrain - dle nastaveni constraint
// interval - constraint - jako v cmp u Qsortu
// -1 vetsi nez interval (min_dt,max_dt)
// +1 mensi
// 0 OK

int TimeGovernor::set_upper_constraint (double upper, std::string message)
{
    if (upper_constraint_ < upper) {
        //do not change upper_constraint_
        return -1;
    } else
        if (lower_constraint_ > upper) {
            // set upper constraint to the lower constraint
            upper_constraint_ = lower_constraint_;
			upper_constraint_message_ = "Forced lower constraint. " + message;
            return 1;
        } else {
            //change upper_constraint_ to upper
            upper_constraint_ = upper;
            upper_constraint_message_ = message;
            return 0;
        }
}



int TimeGovernor::set_lower_constraint (double lower, std::string message)
{   
    if (upper_constraint_ < lower) {
        // set lower constraint to the upper constraint
        lower_constraint_ = upper_constraint_;
        return -1;
    } else
        if (lower_constraint_ > lower) {
            //do not change lower_constraint_
            return 1;
        } else {
            //change lower_constraint_ to lower
            lower_constraint_ = lower;
	        lower_constraint_message_ = message;
            return 0;
        }
}




double TimeGovernor::fix_dt_until_mark() {
    if (steady_) return 0.0;
    end_of_fixed_dt_interval_=-inf_time; // release previous fixed interval
    fixed_time_step_ = estimate_dt();
    is_time_step_fixed_ = true;    //flag means fixed step has been set since now
    return end_of_fixed_dt_interval_ = time_marks_.next(*this, equation_fixed_mark_type())->time();
}




void TimeGovernor::add_time_marks_grid(double step, TimeMark::Type mark_type) const
{
	if (end_time() == inf_time) {
		THROW(ExcTimeGovernorMessage()
				<< EI_Message("Missing end time for making output grid required by key 'time_step' of the output stream.\n")
			);
	}

	marks().add_time_marks(init_time_, step, end_time(), mark_type | eq_mark_type_);
	// always add start time and end time
	marks().add(TimeMark(init_time_, mark_type | eq_mark_type_));
	marks().add(TimeMark(end_time(), mark_type | eq_mark_type_));
}


bool TimeGovernor::is_current(const TimeMark::Type &mask) const
{
    TimeMark::Type type = equation_mark_type() | mask;
    return time_marks_.current(step(), type) != time_marks_.end(type);
}



double TimeGovernor::estimate_dt() const {
    if (is_end()) return 0.0;
    
    if (this->step().lt(end_of_fixed_dt_interval_))    return fixed_time_step_;

    // jump to the first future fix time
    TimeMarks::iterator fix_time_it = time_marks_.next(*this, equation_fixed_mark_type());
    // compute step to next fix time and apply constraints
    double full_step = fix_time_it->time() - t();

    double step_estimate = min(full_step, upper_constraint_);
    step_estimate = max(step_estimate, lower_constraint_); //these two must be in this order

    if (step_estimate == inf_time) return step_estimate;

    // round the time step to have integer number of steps till next fix time
    // this always selects shorter time step,
    // but allows time step larger then constraint by a number close to machine epsilon
    //
    double n_steps = ceil( full_step / step_estimate - time_step_precision);
    step_estimate = full_step / n_steps;
    
    // try to avoid time_step changes
    if (n_steps > 1 && abs(step_estimate - step().length()) < time_step_precision) {
        step_estimate=step().length();
    }
    
    // if the step estimate gets by rounding lower than lower constraint program will not crash
    // will just output a user warning.
    if (step_estimate < lower_constraint_) {
    	static char buffer[1024];
    	sprintf(buffer, "Time step estimate is below the lower constraint of time step. The difference is: %.16f.\n",
                lower_constraint_ - step_estimate);
    	WarningOut() << buffer;
    }
    
    return step_estimate;
}



void TimeGovernor::next_time()
{
    OLD_ASSERT_LE(0.0, t());
    if (is_end()) return;
    

    if (this->step().lt(end_of_fixed_dt_interval_)) {
        // this is done for fixed step
        // make tiny correction of time step in order to avoid big rounding errors
        // tiny correction means that dt_changed 'is NOT changed'
    	if (end_of_fixed_dt_interval_ < inf_time) {
    		fixed_time_step_ = (end_of_fixed_dt_interval_-t()) / round( (end_of_fixed_dt_interval_-t()) / fixed_time_step_ );
    	}

    	recent_steps_.push_front(recent_steps_.front().make_next(fixed_time_step_));


        //checking whether fixed time step has been changed (by fix_dt_until_mark() method) since last time
        if (is_time_step_fixed_)
        {
            is_time_step_fixed_ = false;       
            
            //is true unless new fixed_dt is not equal previous time_step
            time_step_changed_ = (step(-2).length() != step().length());
        }
        else
            time_step_changed_ = false;
    }
    else
    {
        // this is done if end_of_fixed_dt_interval is not set (means it is equal to -infinity)
        double dt=estimate_dt();
        TimeStep step_ = recent_steps_.front().make_next(dt);
        //DebugOut().fmt("step: {}, end: {}\n", step_.length(), step_.end());
        recent_steps_.push_front(step_);
        //DebugOut().fmt("last: {}, new: {}\n",step(-1).length(),step().length());
        time_step_changed_= (step(-2).length() != step().length());
    }

    last_lower_constraint_ = lower_constraint_;
    last_upper_constraint_ = upper_constraint_;
    // refreshing the upper_constraint_
    upper_constraint_ = min(end_time_ - t(), max_time_step_);
    lower_constraint_ = min_time_step_;
    lower_constraint_message_ = "Permanent minimal constraint, in next time.";
    upper_constraint_message_ = "Permanent maximal constraint, in next time.";
}


double TimeGovernor::reduce_timestep(double factor) {
    double prior_dt = dt();
    double new_upper_constraint = factor * dt();

    // Revert time.
//    DebugOut().fmt("tg idx: {}\n", recent_steps_.front().index());
    recent_steps_.pop_front();
//    DebugOut().fmt("tg idx: {}\n", recent_steps_.back().index());
    upper_constraint_ = last_upper_constraint_;
    lower_constraint_ = last_lower_constraint_;

    // Set constraint.
    int current_minus_new = set_upper_constraint(new_upper_constraint, "Reduce time step.");
    if (current_minus_new < 0)
        // current constraint < reduced dt, should not happen
        THROW(ExcMessage() << EI_Message("Internal error."));

    next_time();

    // Return false if we hit lower time step constraint.
    return dt() / prior_dt;
}



const TimeStep &TimeGovernor::step(int index) const {
    unsigned int back_idx;
    if (index < 0) {
        back_idx = static_cast<unsigned int>(-index-1);
    } else {
        back_idx = static_cast<unsigned int>(recent_steps_[0].index() - index);
    }
    if ( back_idx >= recent_steps_.size())
        THROW(ExcMissingTimeStep() << EI_Index(index) << EI_BackIndex(back_idx) << EI_HistorySize(recent_steps_.size()));

    return recent_steps_[back_idx];
}




void TimeGovernor::view(const char *name) const
{
	static char buffer[1024];
#ifdef FLOW123D_DEBUG_MESSAGES
	sprintf(buffer, "TG[%s]:%06d    t:%10.4f    dt:%10.6f    dt_int<%10.6f,%10.6f>    end_time: %f end_fixed_time: %f type: 0x%x\n",
	            name, tlevel(), t(), dt(), lower_constraint_, upper_constraint_, end_time_,  end_of_fixed_dt_interval_, eq_mark_type_);
#else
	sprintf(buffer, "TG[%s]:%06d    t:%10.4f    dt:%10.6f    dt_int<%10.6f,%10.6f>\n",
	            name, tlevel(), t(), dt(), lower_constraint_, upper_constraint_ );
#endif
	MessageOut() << buffer;
    MessageOut().fmt("Lower time step constraint [{}]: {} \nUpper time step constraint [{}]: {} \n",
            lower_constraint_, lower_constraint_message_.c_str(), 
            upper_constraint_, upper_constraint_message_.c_str() );
}




ostream& operator<<(ostream& out, const TimeGovernor& tg)
{
	static char buffer[1024];
	sprintf(buffer, "\n%06d    t:%10.4f    dt:%10.6f    dt_int<%10.6f,%10.6f>\n",
            tg.tlevel(), tg.t(), tg.dt(), tg.lower_constraint(), tg.upper_constraint());
	return (out << buffer);
}


