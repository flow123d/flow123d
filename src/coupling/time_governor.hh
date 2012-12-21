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

#include "system/global_defs.h"
#include "system/system.hh"
#include "time_marks.hh"

namespace Input {
    class Record;
    namespace Type {
        class Record;
    }
}

/**
 * @brief
 * Basic time management functionality for unsteady (and steady) solvers (class Equation).
 *
 * This class provides algorithm for selecting next time step, and information about current time step frame.
 * Step estimating is constrained by several bounds (permanent maximal and minimal time step, upper 
 * and lower constraint of time step). The permanent constraints are set in the constructor from the input 
 * record so that user can set the time step constraints for the whole model. 
 * Function set_permanent_constraint() should be used only in very specific cases and possibly right after 
 * the constructor before using other functions of TimeGovernor.
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
 * Information provided by time governor includes:
 * - actual time, last time, end time
 * - actual time step
 * - number of the time level
 * - end of interval with fixed time step
 * - time comparison
 * - static pointer to time marks
 *
 * TODO: better implementation of steady TimeGovernor, be careful with infinite values namely where we do some calculations with time:
 * - TimeMarks::is_current
 * - TimeGovernor:: lt le ge gt
 *
 * TODO:
 * - still we have problems with time comparisons
 * 1) TimeMarks can merge marks only with fixed precision, since they are shared by several equations with possibly different timesteps
 * 2) queries
 *
 *
 *
 */

class TimeGovernor
{
public:
    
    /**
     * @brief Constructor for unsteady solvers.
     *
     * @param input accessor to input data
     * @param fixed_time_mask TimeMark mask used to select fixed time marks from all the time marks. 
     * This value is bitwise added to the default one defined in TimeMarks::type_fixed_time().
     *
     */
   TimeGovernor(const Input::Record &input,
                const TimeMark::Type fixed_time_mask = 0x0);

   /**
    * @brief Constructor - steady time governor.
    *
    * Optionally you can set initial time for "one step" steady problems.
    * @see default constructor.
    * 
    */
   TimeGovernor(double init_time);

   /**
    * @brief Deafult constructor - steady time governor.
    * 
    * We can have "zero step" steady problem (no computation, e.g. EquationNothing) and one step steady problem
    * (e.g. steady water flow).
    * 
    * Time is set to zero, time step and end time to infinity.
    * 
    * First call of next_time() push the actual time to infinity.
    * 
    * However, you have to use full constructor for the "steady problem" that have time variable input data.
    * 
    * Has a private pointer to static TimeMarks and can access them by marks().
    */
   TimeGovernor();

   static Input::Type::Record input_type;

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
     * @return -1, 0 or 1 according to the succes.
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
    inline double fix_dt_until_mark() {
        if (steady) return inf_time;
        end_of_fixed_dt_interval=-inf_time; // release previous fixed interval
        fixed_dt = estimate_dt();
        dt_fixed_now = true;    //flag means fixed step has been set since now
        return end_of_fixed_dt_interval = time_marks->next(*this, fixed_time_mark_mask)->time();
    }

    /**
     * @brief Proceed to the next time according to current estimated time step.
     */
    void next_time();

    /**
     * Getter for time marks.
     */
    static inline TimeMarks &marks()
            {return *time_marks;}

    /**
     * Simpler interface to TimeMarks::is_current().
     */
    inline bool is_current(const TimeMark::Type &mask) const
        {return time_marks->is_current(*this, mask); }

    /**
     * Simpler interface to TimeMarks::next().
     */
    inline TimeMarks::iterator next(const TimeMark::Type &mask) const
        {return time_marks->next(*this, mask);}

    /**
     * Simpler interface to TimeMarks::last().
     */
    inline TimeMarks::iterator last(const TimeMark::Type &mask) const
        {return time_marks->last(*this, mask);}

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
        {return end_of_fixed_dt_interval;}

    /**
     *  Getter for dt_changed. Returns whether the time step has been changed.
     */
    inline bool is_changed_dt() const
        {return dt_changed;}


    /**
     * End of actual time interval; i.e. where the solution is computed.
     */
    inline double t() const
        {return time;}

    /**
     * Previous time step.
     */
    inline double last_dt() const
        {return last_time_step;}

    /**
     * Length of actual time interval; i.e. the actual time step.
     */
    inline double dt() const
        {return time_step;}

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
        {return time+estimate_dt();}

    /// End time.
    inline double end_time() const
    { return end_time_; }

    /// Returns true if the actual time is greater than or equal to the end time.
    inline bool is_end() const
        { return (this->ge(end_time_) || time == inf_time); }
        
    /// Returns true if the time governor is used for steady problem.
    inline double is_steady() const
    { return steady; }

    /**
     * Performs rounding safe comparison time > other_time, i.e. time is strictly greater than given parameter
     * other_time with precision relative to the magnitude of the numbers time step.
     * TODO: introduce type TimeDouble with overloaded comparison operators, use it consistently in TimeMarks.
     */
    inline bool gt(double other_time) const
        {
            return ! (time <= other_time
            + 16*numeric_limits<double>::epsilon()*max(abs(time),abs(other_time)) );
        }

    /**
     * Performs rounding safe comparison time >= other_time See @fn gt
     */
    inline bool ge(double other_time) const
    {
        return time >= other_time
        - 16*numeric_limits<double>::epsilon()*max(abs(time),abs(other_time));
    }

    /**
     * Performs rounding safe comparison time < other_time. See @fn gt
     */
    inline bool lt(double other_time) const
    {
        double b=other_time
                - 16*numeric_limits<double>::epsilon()*max(abs(time),abs(other_time));
        //DBGMSG("time: %e otime: %e eps: %e result: %d\n", time, b,
        //        time - b,
        //         time >= b);
        return ! (time >= b);
    }

    /**
     * Performs rounding safe comparison time <= other_time. See @fn gt
     */
    inline bool le(double other_time) const
    {
        return time <= other_time
        + 16*numeric_limits<double>::epsilon()*max(abs(time),abs(other_time));
    }

    /**
     * Returns the time level.
     */
    inline int tlevel() const
        {return time_level;}

    /**
     * Prints out TimeGovernor status -- time level, end time, actual time and step.
     */
    void view() const
    {
        xprintf(MsgDbg, "\nTG: level: %d end_time: %f time: %f step: %f upper: %f lower: %f end_fixed_time: %f\n",time_level, end_time_, time, time_step, upper_constraint_, lower_constraint_, end_of_fixed_dt_interval);
        //xprintf(Msg, "TG: level: %d end_time: %f time: %f step: %f upper: %f lower: %f end_fixed_time: %f\n",time_level, end_time_, time, time_step, upper_constraint_, lower_constraint_, end_of_fixed_dt_interval);
    }

    /// Infinity time used for steady case.
    static const double inf_time;

private:
    inline double comparison_fracture() const
    {
        if (time_level!=0 && time_step <=numeric_limits<double>::max() ) return comparison_precision * time_step;
        else return numeric_limits<double>::epsilon();
    }


    /// We consider time difference is zero if it is less then comparison_precision * time_step.
    static const double comparison_precision;
    /// Technical bound for the time step given by finite precision.
    static const double time_step_lower_bound;
    /// Rounding precision for computing number of steps. Used in estimate_dt().
    static const double round_n_steps_precision;

    /// Number of time_next calls, i.e. total number of performed time steps.
    int time_level;
    /// End of actual time interval; i.e. where the solution is computed.
    double time;
    /// Beginning of the actual time interval; i.e. the time of last computed solution.
    //double last_time;
    /// End of interval if fixed time step.
    double end_of_fixed_dt_interval;
    /// End time of the simulation.
    double end_time_;

    /// Length of actual time interval; i.e. the actual time step.
    double time_step;
    /// Time step just before last_time.
    double last_time_step;
    /// Next fixed time step.
    double fixed_dt;
    /// Flag that is set when the fixed step is set (lasts only one time step).
    bool dt_fixed_now;
    /// Flag is set if the time step has been changed (lasts only one time step).
    bool dt_changed;

    
    /// Upper constraint for the choice of the next time step. Relaxed after every dt choice. OBSOLETE
    //double time_step_constraint;
    
    /// Upper constraint for the choice of the next time step.
    double upper_constraint_;
    /// Lower constraint for the choice of the next time step.
    double lower_constraint_;
    /// Permanent upper limit for the time step.
    double max_time_step;
    /// Permanent lower limit for the time step.
    double min_time_step;

    /**
     * When the next time is chosen we need only the lowest fix time. Therefore we use
     * minimum priority queue of doubles based on the vector container.
     * This is one global set of time marks for the whole problem and is shared umong all equations. 
     * Therefore this object is static constant pointer.
     */
    static TimeMarks * const time_marks;
    
    /// TimeMark type that masks the fixed time mark. It is set by constructor (unsteady case).
    const TimeMark::Type fixed_time_mark_mask;
    
    /// True if the time governor is used for steady problem.
    const bool steady;

};

#endif /* TIME_HH_ */
