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
#include "unit_si.hh"

/*******************************************************************
 * implementation of TimeGovernor static values and methods
 */

//initialize constant pointer to TimeMarks object
TimeMarks TimeGovernor::time_marks_ = TimeMarks();

//const double TimeGovernor::time_step_lower_bound = numeric_limits<double>::epsilon();

#define MAX_END_TIME 5.0e+17
#define MAX_END_TIME_STR "5.0e+17"

const double TimeGovernor::inf_time =  numeric_limits<double>::infinity();
const double TimeGovernor::max_end_time = MAX_END_TIME; // more then age of universe in seconds.
const double TimeGovernor::time_step_precision = 16*numeric_limits<double>::epsilon();


using namespace Input::Type;


const Tuple & TimeGovernor::get_input_time_type(double lower_bound, double upper_bound)
{
    return Tuple("TimeValue", "A time with optional unit specification.")
        .declare_key("time", Double(lower_bound, upper_bound), Default::obligatory(),
                                    "The time value." )
		.declare_key("unit", String(), Default::read_time("Common time unit of the equation's Time Governor."),
				"The time unit. Possible values: 's' seconds, 'min' minutes, 'h' hours, 'd' days, 'y' years.")
		.close();
}


const Record & TimeGovernor::get_input_type() {
    static const Tuple &dt_step =
        Tuple("DtLimits", "Time dependent changes in min_dt and max_dt limits.")
            .declare_key("time", TimeGovernor::get_input_time_type(), Default::obligatory(),
                    "The start time of dt step set.")
            .declare_key("min_dt", TimeGovernor::get_input_time_type(), Default::read_time("'min_dt' value of TimeGovernor."),
                    "Soft lower limit for the time step.")
            .declare_key("max_dt", TimeGovernor::get_input_time_type(), Default::read_time("'max_dt' value of TimeGovernor."),
                    "Whole time of the simulation if specified, infinity else.")
            .close();

    return Record("TimeGovernor",
            "Time axis settings of the simulation.\n"
            "The settings is specific to a particular equation.\n"
    		"TimeGovernor allows to:\n"
    		" - define start time and end time of simulation\n"
    		" - define lower and upper limits of time steps\n"
    		" - direct fixed time marks of whole simulation\n"
    		" - set global time unit of equation (see 'common_time_unit' key)\n"
    		"Limits of time steps are defined by keys 'min_dt', 'max_dt', 'init_dt' and 'dt_limits'. Key "
    		"'init_dt' has the highest priority and allows set fix size of time steps. Pair of keys 'min_dt' "
    		"and 'max_dt' define interval of time steps. Both previous cases ('init_dt' or pair 'min_dt' "
    		"and 'max_dt') set global limits of whole simulation. In contrasts, 'dt_limits' allow set "
    		"time-dependent function of min_dt/max_dt. Used time steps of simulation can be printed to YAML "
    		"output file (see 'write_used_timesteps'.\n"
    		"Fixed time marks define exact values of time steps. They are defined in:\n"
    		" - start time and end time of simulation\n"
    		" - output times printed to output mesh file\n"
    		" - times defined in 'dt_limits' table (optional, see 'add_dt_limits_time_marks' key)")
		.allow_auto_conversion("max_dt")
		.declare_key("start_time", TimeGovernor::get_input_time_type(), Default("0.0"),
					"Start time of the simulation.")
		.declare_key("end_time", TimeGovernor::get_input_time_type(), Default(MAX_END_TIME_STR),
				"End time of the simulation.\n"
				"The default value is higher than the age of the Universe (given in seconds).")
		.declare_key("init_dt", TimeGovernor::get_input_time_type(0.0), Default("0.0"),
				"Initial guess for the time step.\n"
				"It applies to equations that use an adaptive time stepping. "
				"If set to 0.0, the time step is determined in fully autonomous "
				"way, assuming the equation supports it.")
		.declare_key("min_dt", TimeGovernor::get_input_time_type(0.0),
				Default::read_time("Machine precision."),
				"Soft lower limit for the time step.\n"
				"Equation using an adaptive time stepping cannot suggest smaller time step. "
				"The actual time step can only decrease below the limit in order to match "
				"the prescribed input or output times.")
		.declare_key("max_dt", TimeGovernor::get_input_time_type(0.0),
				Default::read_time("Whole time of the simulation if specified, infinity else."),
				"Hard upper limit for the time step.\n"
				"The actual time step can only increase above the limit in order to match "
				"the prescribed input or output times.")
		.declare_key("dt_limits", Array(dt_step), Default::optional(),
				"Allow to set a time dependent changes in ``min_dt`` and ``max_dt`` limits. This list is processed "
				"at individual times overwriting previous values of ``min_dt``/``max_dt``. Limits equal to 0 are "
				"ignored and replaced with ``min_dt``/``max_dt`` values.")
		.declare_key("add_dt_limits_time_marks", Bool(), Default("false"), "Add all times defined in ``dt_limits`` "
			    "table to the list of fixed TimeMarks.")
		.declare_key("write_used_timesteps", FileName::output(), Default::optional(),
				"Write used time steps to the given file in YAML format corresponding with the format of ``dt_limits``.")
		.declare_key("common_time_unit", String(), Default("\"s\""),
				"Common time unit of the equation.\nThis unit will be used for all time inputs and outputs "
				"within the equation. Individually, the common time unit can be overwritten for every declared time.\n"
				"Time units are used in the following cases:\n"
				"1) Time units of time value keys in: TimeGovernor, FieldDescriptors.\n"
				"   The common time unit can be overwritten for every declared time.\n"
				"2) Time units in: \n"
				"   a) input fields: FieldFE and FieldTimeFunction\n"
				"   b) time steps definition of OutputTimeSet\n"
				"   Common time unit can be overwritten by one unit value for every whole mesh data file or time function.\n"
				"3) Time units in output files: observation times, balance times, frame times of VTK and GMSH\n"
				"   Common time unit cannot be overwritten in these cases."
				)
		.close();
}



/*******************************************************************
 * implementation of TimeUnitConversion
 */

TimeUnitConversion::TimeUnitConversion(std::string user_defined_unit)
: unit_string_(user_defined_unit)
{
    coef_ = UnitSI().s().convert_unit_from(user_defined_unit);
}



TimeUnitConversion::TimeUnitConversion()
: coef_(1.0), unit_string_("s") {}



double TimeUnitConversion::read_time(Input::Iterator<Input::Tuple> time_it, double default_time) const {
	if (time_it) {
	    double time = time_it->val<double>("time");
	    string time_unit;
		if (time_it->opt_val<string>("unit", time_unit)) {
			return ( time * UnitSI().s().convert_unit_from(time_unit) );
		} else {
			return ( time * coef_ );
		}
	} else {
		ASSERT(default_time!=std::numeric_limits<double>::quiet_NaN()).error("Undefined default time!");
		return default_time;
	}
}



double TimeUnitConversion::read_coef(Input::Iterator<string> unit_it) const {
	if (unit_it) {
		return UnitSI().s().convert_unit_from(*unit_it);
	} else {
		return coef_;
	}
}



/*******************************************************************
 * implementation of TimeStep
 */

TimeStep::TimeStep(double init_time, std::shared_ptr<TimeUnitConversion> time_unit_conversion) :
index_(0),
length_(1.0),
end_(init_time),
time_unit_conversion_(time_unit_conversion)
{}



TimeStep::TimeStep() :
index_(0),
length_(TimeGovernor::inf_time),
end_(-TimeGovernor::inf_time)
{
	time_unit_conversion_ = std::make_shared<TimeUnitConversion>();
}




TimeStep::TimeStep(const TimeStep &other):
index_(other.index_),
length_(other.length_),
end_(other.end_),
time_unit_conversion_(other.time_unit_conversion_)
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
    ts.time_unit_conversion_=time_unit_conversion_;
    return ts;
}



bool TimeStep::safe_compare(double t1, double t0) const
{
    return t1 >= t0
            - 16*numeric_limits<double>::epsilon()*(1.0+max(abs(t1),abs(t0)));
}



double TimeStep::read_time(Input::Iterator<Input::Tuple> time_it, double default_time) const {
	return time_unit_conversion_->read_time(time_it, default_time);
}



double TimeStep::read_coef(Input::Iterator<string> unit_it) const {
	return time_unit_conversion_->read_coef(unit_it);
}



double TimeStep::get_coef() const {
	return time_unit_conversion_->get_coef();
}



ostream& operator<<(ostream& out, const TimeStep& t_step) {
    out << "time: " << t_step.end() << "step: " << t_step.length() << endl;
    return out;
}



/*******************************************************************
 * implementation of TimeGovernor
 */

TimeGovernor::TimeGovernor(const Input::Record &input, TimeMark::Type eq_mark_type, bool timestep_output)
: timestep_output_(timestep_output)
{
    // use new mark type as default
    if (eq_mark_type == TimeMark::none_type) eq_mark_type = marks().new_mark_type();

    try {

        string common_unit_string=input.val<string>("common_time_unit");
        time_unit_conversion_ = std::make_shared<TimeUnitConversion>(common_unit_string);
        limits_time_marks_ = input.val<bool>("add_dt_limits_time_marks");

        // Get rid of rounding errors.
        double end_time = read_time( input.find<Input::Tuple>("end_time") );
        if (end_time> 0.99*max_end_time) end_time = max_end_time;

        // set permanent limits
    	init_common(read_time( input.find<Input::Tuple>("start_time") ),
    				end_time,
    				eq_mark_type);
    	Input::Array limits_array = Input::Array();
    	input.opt_val("dt_limits", limits_array);
		set_dt_limits(
			read_time( input.find<Input::Tuple>("min_dt"), min_time_step_),
			read_time( input.find<Input::Tuple>("max_dt"), max_time_step_),
			limits_array
			);

    	// check key write_used_timesteps, open YAML file, print first time step
        if (timestep_output_)
            if (input.opt_val("write_used_timesteps", timesteps_output_file_) ) {
                try {
                    timesteps_output_file_.open_stream(timesteps_output_);
                } INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, input)
                timesteps_output_ << "- [ " << t() << ", " << dt_limits_table_[0].min_dt << ", " << dt_limits_table_[0].max_dt << " ]\n";
                last_printed_timestep_ = t();
            }

        double init_dt=read_time( input.find<Input::Tuple>("init_dt") );
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
: dt_limits_pos_(0), timestep_output_(false)
{
	time_unit_conversion_ = std::make_shared<TimeUnitConversion>();
	init_common( init_time, inf_time, TimeMark::every_type);
    // fixed time step
    if (dt < time_step_precision)
        THROW(ExcTimeGovernorMessage() << EI_Message("Fixed time step smaller then machine precision. \n") );

    fixed_time_step_=dt;
    is_time_step_fixed_=true;
    time_step_changed_=true;
    end_of_fixed_dt_interval_ = inf_time;

    // fill table limits with two records (start time, end time)
    dt_limits_table_.push_back( DtLimitRow(init_time, dt, dt) );
    dt_limits_table_.push_back( DtLimitRow(inf_time, dt, dt) );

    lower_constraint_=min_time_step_=dt;
    lower_constraint_message_ = "Initial time step set by user.";
    upper_constraint_=max_time_step_=dt;
    upper_constraint_message_ = "Initial time step set by user.";
    
    //time_step_=dt;
}


// steady time governor constructor
TimeGovernor::TimeGovernor(double init_time, TimeMark::Type eq_mark_type)
: dt_limits_pos_(0), timestep_output_(false)
{
    // use new mark type as default
    if (eq_mark_type == TimeMark::none_type) eq_mark_type = marks().new_mark_type();

    time_unit_conversion_ = std::make_shared<TimeUnitConversion>();
	init_common(init_time, inf_time, eq_mark_type);

	// fill table limits with two records (start time, end time)
    dt_limits_table_.push_back( DtLimitRow(init_time, min_time_step_, max_time_step_) );
    dt_limits_table_.push_back( DtLimitRow(inf_time, min_time_step_, max_time_step_) );

	steady_ = true;
}


TimeGovernor::~TimeGovernor()
{
	if ( !(timesteps_output_file_ == FilePath()) && timestep_output_ ) {
		timesteps_output_.close();
	}
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
    recent_steps_.push_front(TimeStep(init_time, time_unit_conversion_));
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
   		max_time_step_=upper_constraint_= end_time - init_time_;
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
    //if (end_time_ != inf_time)
    time_marks_.add( TimeMark(end_time_, equation_fixed_mark_type()) );
}



void TimeGovernor::set_dt_limits( double min_dt, double max_dt, Input::Array dt_limits_list)
{
	dt_limits_table_.clear();

	// check min_dt and max_dt set by user
	if (min_dt < time_step_precision) {
		THROW(ExcTimeGovernorMessage() << EI_Message("'min_dt' smaller than machine precision.\n") );
	}
	if (max_dt < min_dt) {
		THROW(ExcTimeGovernorMessage() << EI_Message("'max_dt' smaller than 'min_dt'.\n") );
	}

	bool first_step = true;
    if (dt_limits_list.size())
    	for(auto it = dt_limits_list.begin<Input::Tuple>(); it != dt_limits_list.end(); ++it) {
			double time = read_time( it->find<Input::Tuple>("time"));

			if (first_step) { // we need special check before setting first time step to the table
				if (time > init_time_) { // table starts later than simulation, we need to add start time to the table
					dt_limits_table_.push_back( DtLimitRow(init_time_, min_dt, max_dt) );
				}
				first_step = false;
			}

			// next cases will be skipped
			if (time < init_time_) {
				WarningOut().fmt("Time {} define in 'dt_limits' table at address {} is lesser than start time of simulation "
						"and can be skipped.\n", time, dt_limits_list.address_string());
			}
			if (dt_limits_table_.size() && (time <= dt_limits_table_[dt_limits_table_.size()-1].time) ) {
				WarningOut().fmt("Time {} define in 'dt_limits' table at address {} is in incorrect order "
						"and will be skipped.\n", time, dt_limits_list.address_string());
				continue;
			}
			if ((time > end_time_) ) {
				WarningOut().fmt("Time {} define in 'dt_limits' table at address {} is greater than end time of simulation "
						"and will be skipped.\n", time, dt_limits_list.address_string());
				continue;
			}

			double min = read_time( it->find<Input::Tuple>("min_dt"), 0.0);
			if (min == 0.0) min = min_dt;
			double max = read_time( it->find<Input::Tuple>("max_dt"), 0.0);
			if (max == 0.0) max = max_dt;

			if (min < time_step_precision) {
				THROW(ExcTimeGovernorMessage() << EI_Message("'min_dt' in 'dt_limits' smaller than machine precision.\n") );
			}
			if (max < min) {
				THROW(ExcTimeGovernorMessage() << EI_Message("'max_dt' in 'dt_limits' smaller than 'min_dt'.\n") );
			}

			dt_limits_table_.push_back( DtLimitRow(time, min, max) );
			if (limits_time_marks_) this->marks().add(TimeMark(time, this->equation_fixed_mark_type()));
		}

    if (dt_limits_table_.size() == 0) {
    	// add start time to limit table if it is empty
    	dt_limits_table_.push_back( DtLimitRow(init_time_, min_dt, max_dt) );
    }
    if (dt_limits_table_[dt_limits_table_.size()-1].time < end_time_) {
    	// add time == end_time_ to limits table, we need only for check time, not for limits
    	dt_limits_table_.push_back( DtLimitRow(end_time_, min_dt, max_dt) );
    }

    dt_limits_pos_ = 0;
    while (dt_limits_table_[dt_limits_pos_+1].time <= init_time_) {
    	++dt_limits_pos_;
    }

    set_permanent_constraint();
}


void TimeGovernor::set_permanent_constraint()
{
    lower_constraint_ = min_time_step_ = max(dt_limits_table_[dt_limits_pos_].min_dt, time_step_precision);
    lower_constraint_message_ = "Permanent minimal constraint, custom.";
    upper_constraint_ = max_time_step_ = min(dt_limits_table_[dt_limits_pos_].max_dt, end_time_-t());
    upper_constraint_message_ = "Permanent maximal constraint, custom.";
    ++dt_limits_pos_;
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
        DebugOut().fmt("Time step estimate is below the lower constraint of time step. The difference is: {:.16f}.\n",
                lower_constraint_ - step_estimate);
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

    if (step().end() >= dt_limits_table_[dt_limits_pos_].time) set_permanent_constraint();

	// write time step to YAML file
    if ( !(timesteps_output_file_ == FilePath()) && timestep_output_ ) {
    	double time = t();
    	if (time > last_printed_timestep_) {
    		if (is_end()) timesteps_output_ << "- [ " << time << ", 0, 0 ]\n";
    		else timesteps_output_ << "- [ " << time << ", " << lower_constraint_ << ", " << upper_constraint_ << " ]\n";
    		last_printed_timestep_ = time;
    	}
	}
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
#ifdef FLOW123D_DEBUG_MESSAGES
    MessageOut().fmt(
            "TG[{}]:{:06d}    t:{:10.4f}    dt:{:10.6f}    dt_int<{:10.6f},{:10.6f}>    "
            "end_time: {} end_fixed_time: {} type: {:#x}\n",
            name, tlevel(), t(), dt(), lower_constraint_, upper_constraint_,
            end_time_,  end_of_fixed_dt_interval_, eq_mark_type_.bitmap_);

    MessageOut().fmt("Lower time step constraint [{}]: {} \nUpper time step constraint [{}]: {} \n",
            lower_constraint_, lower_constraint_message_.c_str(), 
            upper_constraint_, upper_constraint_message_.c_str() );
#else
    MessageOut().fmt(
            "TG[{}]:{:06d}    t:{:10.4f}    dt:{:10.6f}    dt_int<{:10.6f},{:10.6f}>\n",
            name, tlevel(), t(), dt(), lower_constraint_, upper_constraint_);
#endif
}



double TimeGovernor::read_time(Input::Iterator<Input::Tuple> time_it, double default_time) const {
	return time_unit_conversion_->read_time(time_it, default_time);
}



double TimeGovernor::read_coef(Input::Iterator<string> unit_it) const {
	return time_unit_conversion_->read_coef(unit_it);
}



double TimeGovernor::get_coef() const {
	return time_unit_conversion_->get_coef();
}



string TimeGovernor::get_unit_string() const {
	return time_unit_conversion_->get_unit_string();
}




ostream& operator<<(ostream& out, const TimeGovernor& tg)
{
	static char buffer[1024];
	sprintf(buffer, "\n%06d    t:%10.4f    dt:%10.6f    dt_int<%10.6f,%10.6f>\n",
            tg.tlevel(), tg.t(), tg.dt(), tg.lower_constraint(), tg.upper_constraint());
	return (out << buffer);
}


