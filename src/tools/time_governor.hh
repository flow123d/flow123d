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
 * $Id: time_governor.hh 1974 2012-11-14 13:46:11Z pavel.exner $
 * $Revision: 1974 $
 * $LastChangedBy: pavel.exner $
 * $LastChangedDate: 2012-11-14 14:46:11 +0100 (St, 14 lis 2012) $
 *
 * @file
 * @brief Basic time management class.
 *  @author Jan Brezina
 */

#ifndef TIME_HH_
#define TIME_HH_

#include <limits>
#include <cmath>
#include <algorithm>
#include <boost/circular_buffer.hpp>


#include "system/global_defs.h"
#include "system/system.hh"
#include "input/accessors.hh"
#include "tools/time_marks.hh"

namespace Input {
    class Record;
    namespace Type {
        class Record;
    }
}




/**
 * Structure to store and provide information about
 * one time step.
 */

class TimeStep {
public:
    /**
     * Constructor of the zero time step.
     */
    TimeStep(double init_time);

    /**
     * Deafult constructor.
     */
    TimeStep();

    /**
     * Copy constructor.
     */
    TimeStep(const TimeStep &other);

    /**
     * Create subsequent time step.
     */
    TimeStep make_next(double new_length) const;

    /**
     * Create subsequent time step, with the @end_time
     * explicitly specified. This allow slight discrepancy to
     * overcome rounding errors in the case of fixed time step.
     * Otherwise using small fixed time step, we may miss long term fixed
     * goal time.
     *
     */
    TimeStep make_next(double new_lenght, double end_time) const;

    /**
     * Getters.
     */
    unsigned int index() const {return index_;}
    double length() const { return length_;}
    double end() const { return end_;}
    double begin() const {return begin_;}

private:
    /// Index of the step is index if the end time. Zero time step is artificial.
    unsigned int index_;
    /// Length of the time step. Theoretically @p end minus @p begin.
    /// However may be slightly different due to rounding errors.
    double length_;
    /// Beginning of the time step (equal to the
    double begin_;
    /// End time point of the time step.
    double end_;
};



/**
 * @brief
 * Basic time management functionality for unsteady (and steady) solvers (class Equation).
 *
 * <h2> Common features and unsteady time governor (TG) </h2>
 * 
 * This class provides algorithm for selecting next time step, and information about current time step frame.
 * Step estimating is constrained by several bounds (permanent maximal and minimal time step, upper 
 * and lower constraint of time step). The permanent constraints are set in the constructor from the input 
 * record so that user can set the time step constraints for the whole model. 
 * Function set_permanent_constraint() should be used only in very specific cases and possibly right after 
 * the constructor before using other functions of TG.
 * 
 * Choice of the very next time step can be constrained using functions set_upper_constraint() and set_lower_constraint(). 
 * Lower and upper constraints are set equal to permanent ones in the constructor and can only 
 * become stricter. If one tries to set these constraints outside the interval of the previous constraints,
 * nothing is changed and a specified value is returned. Upper and lower constraints are reset in function
 * next_time() to the permanent constraints.
 * 
 * 
 * The later one can be called multiple times with various
 * constraint values and we use the minimum of them. Function next_time() choose the next time step in such a way that it meets actual constraints and
 * a uniform discrete time grid with this step hits the nearest fixed time in lowest possible number of steps.
 *
 * The fixed times are time marks of TimeMarks object passed at construction time with particular mask.
 *
 * There is just one set of time marks for the whole problem. Therefore TimeMarks object is static and is shared umong
 * all the equations and time governors. Each equation creates its own specific time mark type.
 * 
 * Solution of the next time step proceeds in following steps:
 * -# Estimate next time step and set constraints (this can be usefull in hc_seq_explicit for estimating next transport time)
 * -# Fix next time step up to the next time mark (this is necessary for ConvectionTransport)
 * -# Proceed to the next time when solution is available. (this can replace solved flag in equation classes)
 *
 * Information provided by TG includes:
 * - actual time, last time, end time
 * - actual time step
 * - number of the time level
 * - end of interval with fixed time step
 * - time comparison
 * - static pointer to time marks
 * 
 * <h2> Steady time governor</h2>
 * 
 * Steady TG can be constructed by default constructor (initial time is zero) or by 
 * constructor with initial time as parameter. End time and time step are set to infinity. 
 * One can check if the time governor is steady by calling is_steady(). 
 * Calling estimate_dt() will return infinity.
 * 
 * Setting constraints have no consequences. Calling fix_dt_until_mark() will only return zero 
 * and will not do anything.
 * 
 * The steady TG works in two states. At first the time is set to initial and time level 
 * is equal zero. To use steady TG properly one should call next_time() after the computation 
 * of steady problem is done. Current time is then set to infinity, time level is set to 1 and 
 * calling estimate_dt() will return zero.
 * 
 * Note: For example class TransportNothing (which computes really nothing) uses also steady TG but
 * it calls next_time() immediately after TG's construction. This means that the 'computation'of transport 
 * is done.
 * 
 *
 * TODO:
 * - still we have problems with time comparisons
 * 1) TimeMarks can merge marks only with fixed precision, since they are shared by several equations with possibly different timesteps.
 * 2) On the other hand comparing of times by time governor can be done relatively to the current time step.
 *
 *
 */

class TimeGovernor
{
public:

	DECLARE_INPUT_EXCEPTION(ExcTimeGovernorMessage, << EI_Message::val);

    static Input::Type::Record input_type;

    /**
     * Getter for time marks.
     */
    static inline TimeMarks &marks()
            {return time_marks_;}
    

    /**
     * @brief Constructor for unsteady solvers.
     *
     * @param input accessor to input data
     * @param fixed_time_mask TimeMark mask used to select fixed time marks from all the time marks. 
     * This value is bitwise added to the default one defined in TimeMarks::type_fixed_time().
     *
     */
   TimeGovernor(const Input::Record &input,
                TimeMark::Type fixed_time_mask = TimeMark::none_type);


   /**
    * @brief Default constructor - steady time governor.
    * 
    * We can have "zero step" steady problem (no computation, e.g. EquationNothing) and one step steady problem
    * (e.g. steady water flow).
    * 
    * Time is set to zero, time step and end time to infinity.
    * 
    * First call of next_time() pushes the actual time to infinity.
    * 
    * However, you have to use full constructor for the "steady problem" that has time-variable input data.
    * 
    * Has a private pointer to static TimeMarks and can access them by marks().
    */
   explicit TimeGovernor(double init_time=0.0,
		   	    TimeMark::Type fixed_time_mask = TimeMark::none_type);

   /**
    * The aim of this constuctor is simple way to make a time governor without Input interface.
    *
    * TODO: Partially tested as part of field test. Needs its own unit test.
    */
   TimeGovernor(double init_time, double dt);




   /**
    * @brief Sets permanent constraints for time step.
    * 
    * This function should not be normally used. These values are to be set in constructor
    * from the input record or by default.
    * @param min_dt is the minimal value allowed for time step
    * @param max_dt is the maximal value allowed for time step
    */
   void set_permanent_constraint( double min_dt, double max_dt);

    
    /**
     * @brief Sets upper constraint for the next time step estimating.
     * 
     * This function can only make the constraint stricter. Upper constraint is reset to @p max_dt in next_time().
     * @param upper is the upper constraint for time step
     * @return -1, 0 or 1 according to the success
     * 
     * - -1: constraint is higher than the permanent upper constraint @p max_dt. Setting failed, no change happened.
     * - 0: constraint is in the interval of permanent constraints @p min_dt and @p max_dt. The upper constraint has been set.
     * - 1: constraint is lower than permanent lower constraint @p min_dt. Setting failed, no change happened.
     */
    int set_upper_constraint(double upper);
    
    /**
     * @brief Sets lower constraint for the next time step estimating. 
     * @return -1, 0 or 1 according to the success.
     * @see set_upper_constrain().
     */
    int set_lower_constraint(double lower);
    
    /**
     * @brief Fixing time step until fixed time mark.
     * 
     * Fix time step until first fixed time mark. When called inside an already fixed interval, 
     * it overwrites previous setting. 
     * @return actual end of fixed time step.
     */
    double fix_dt_until_mark();

    /**
     * @brief Proceed to the next time according to current estimated time step.
     */
    void next_time();


    /**
     *  Returns reference to required time step in the
     *  recent history. Without parameter the actual time step is returned.
     *  To get previous time steps you have to use negative values of  @p index.
     *  This is to provide better readability.
     */
    inline const TimeStep &step(int index=0) const {
        ASSERT_LE(index, 0);
        unsigned int back_idx = static_cast<unsigned int>(-index);
        ASSERT_LESS(back_idx, recent_steps_.size());
        return recent_steps_[back_idx];
    }

    /**
     *	Specific time mark of the equation owning the time governor.
     */
    inline TimeMark::Type equation_mark_type() const
    { return eq_mark_type_;}

    /**
     *	Specific time mark of the fixed times of the equation owning the time governor.
     */
    inline TimeMark::Type equation_fixed_mark_type() const
    { return eq_mark_type_ | marks().type_fixed_time(); }

    /**
     * Add sequence of time marks starting from the initial time up to the end time with given @p step.
     * Time marks type combines given mark_type (none by default) and native mark type of the time governor.
     */
    void add_time_marks_grid(double step, TimeMark::Type mark_type= TimeMark::none_type) const;

    /**
     * Simpler interface to TimeMarks::is_current().
     */
    inline bool is_current(const TimeMark::Type &mask) const
        {return time_marks_.is_current(*this, equation_mark_type() | mask); }

    /**
     * Simpler interface to TimeMarks::next().
     */
    inline TimeMarks::iterator next(const TimeMark::Type &mask) const
        {return time_marks_.next(*this, mask);}

    /**
     * Simpler interface to TimeMarks::last().
     */
    inline TimeMarks::iterator last(const TimeMark::Type &mask) const
        {return time_marks_.last(*this, mask);}

    /**
     *  Getter for upper constrain.
     */
    inline double upper_constraint() const
        {return upper_constraint_;}
    
    /**
     *  Returns lower constraint.
     */
    inline double lower_constraint() const
        {return lower_constraint_;}
        
    /**
     * End of interval with currently fixed time step. Can be changed by next call of method fix_dt_until_mark.
     */
    inline double end_of_fixed_dt() const
        {return end_of_fixed_dt_interval_;}

    /**
     *  Getter for dt_changed. Returns whether the time step has been changed.
     */
    inline bool is_changed_dt() const
        {return time_step_changed_;}


    /**
     * End of actual time interval; i.e. where the solution is computed.
     */
    inline double t() const
        {return step().end();}

    /**
     * Previous time step.
     */
    inline double last_dt() const
        {if (step().index() >0) return step(-1).length();
        else return inf_time;
        }

    /**
     * Previous time.
     */
    inline double last_t() const
        { return step().begin(); }


    /**
     * Length of actual time interval; i.e. the actual time step.
     */
    inline double dt() const
        {return step().length();}

    /**
     * @brief Estimate choice of next time step according to actual setting of constraints.
     * 
     * Precedence of constraints:
     *
     *  -# meet next fixed time (always satisfied)
     *  -# permanent upper constraint (always satisfied)
     *  -# upper constraint (always satisfied)
     *  -# lower constraint (satisfied if not in conflict with 1.)
     *  -# permanent lower constraint (satisfied if 4.)
     *  -# else writes the difference between lower constraint and estimated time step
     */
    double estimate_dt() const;

    /**
     * Estimate next time.
     */
    inline double estimate_time() const
        {return t()+estimate_dt();}

    /// End time.
    inline double end_time() const
    { return end_time_; }

    /// Returns true if the actual time is greater than or equal to the end time.
    inline bool is_end() const
        { return (this->ge(end_time_) || t() == inf_time); }
        
    /// Returns true if the time governor is used for steady problem.
    inline bool is_steady() const
    { return steady_; }

    /**
     * Performs rounding safe comparison time > other_time, i.e. time is strictly greater than given parameter
     * other_time with precision relative to the magnitude of the numbers time step.
     * TODO: introduce type TimeDouble with overloaded comparison operators, use it consistently in TimeMarks.
     */
    inline bool gt(double other_time) const
        {
            return ! (t() <= other_time
            + 16*numeric_limits<double>::epsilon()*(1.0+max(abs(t()),abs(other_time))) );
        }

    /**
     * Performs rounding safe comparison time >= other_time See @fn gt
     */
    inline bool ge(double other_time) const
    {
        return t() >= other_time
        - 16*numeric_limits<double>::epsilon()*(1.0+max(abs(t()),abs(other_time)));
    }

    /**
     * Performs rounding safe comparison time < other_time. See @fn gt
     */
    inline bool lt(double other_time) const
    {
        double b=other_time
                - 16*numeric_limits<double>::epsilon()*(1.0+max(abs(t()),abs(other_time)));
        return ! (t() >= b);
    }

    /**
     * Performs rounding safe comparison time <= other_time. See @fn gt
     */
    inline bool le(double other_time) const
    {
        return t() <= other_time
        + 16*numeric_limits<double>::epsilon()*(1.0+max(abs(t()),abs(other_time)));
    }

    /**
     * Returns the time level.
     */
    inline int tlevel() const
        {return step().index();}

    /**
     * Prints output of TimeGovernor.
     * @param name is the name of time governor that you want to show up in output (just for your convenience)
     *
     */
    void view(const char *name="") const;

    /// Infinity time used for steady case.
    static const double inf_time;

private:

    /**
     * \brief Common part of the constructors. Set most important parameters, check they are valid and set default values to other.
     *
     * Set main parameters to given values.
     * Check they are correct.
     * Set soft and permanent constrains to the same, the least restricting values.
     * Set time marks for the start time and end time (if finite).
     */
    void init_common(double init_time, double end_time, TimeMark::Type type);


    /// Technical bound for the time step given by finite precision.
    static const double time_step_lower_bound;
    /// Rounding precision for computing number of steps. Used in estimate_dt().
    static const double round_n_steps_precision;

    /// Number of time_next calls, i.e. total number of performed time steps.
    //int time_level_;

    ///
    boost::circular_buffer<TimeStep> recent_steps_;
    /// Initial time.
    double init_time_;
    /// End of actual time interval; i.e. where the solution is computed.
    //double time_;
    /// Beginning of the actual time interval; i.e. the time of last computed solution.
    //double last_time_;
    /// End of interval if fixed time step.
    double end_of_fixed_dt_interval_;
    /// End time of the simulation.
    double end_time_;

    /// Length of actual time interval; i.e. the actual time step.
    //double time_step_;
    /// Time step just before last_time.
    //double last_time_step_;
    /// Next fixed time step.
    double fixed_time_step_;
    /// Flag that is set when the fixed step is set (lasts only one time step).
    bool is_time_step_fixed_;
    /// Flag is set if the time step has been changed (lasts only one time step).
    bool time_step_changed_;

    
    /// Upper constraint for the choice of the next time step.
    double upper_constraint_;
    /// Lower constraint for the choice of the next time step.
    double lower_constraint_;
    /// Permanent upper limit for the time step.
    double max_time_step_;
    /// Permanent lower limit for the time step.
    double min_time_step_;

    /**
     * When the next time is chosen we need only the lowest fix time. Therefore we use
     * minimum priority queue of doubles based on the vector container.
     * This is one global set of time marks for the whole problem and is shared among all equations.
     * Therefore this object is static constant pointer.
     */
    static TimeMarks time_marks_;
    
    /// TimeMark type of the equation.
    TimeMark::Type eq_mark_type_;
    
    /// True if the time governor is used for steady problem.
    bool steady_;

};

/**
 * \brief Redirection operator for TimeGovernor.
 *
 * Currently for debugging purposes.
 * In the future it should be customized for use in combination with
 * streams for various log targets.
 *
 */
ostream& operator<<(ostream& out, const TimeGovernor& tg);


#endif /* TIME_HH_ */
