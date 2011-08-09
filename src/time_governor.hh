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
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Basic time management class.
 */

#ifndef TIME_HH_
#define TIME_HH_

#include <limits>

#include "system/system.hh"
#include "time_marks.hh"

/**
 * @brief
 * Basic time management functionality for unsteady (and steady) solvers (class Equation).
 *
 * This class provides algorithms for selecting next time step, and information about current time step frame.
 * Choice of the next time step can be permanently constrained through function set_permanent_constrain() or one can set
 * constrain only for the very next time step choice through function set_constrain(). The later one can be called multiple times with various
 * constrain values and we use the minimum of them. Function next_time() choose the next time step in such a way that it meets actual constrains and
 * a uniform discrete time grid with this step hits the nearest fixed time in lowest possible number of steps.
 *
 * The fixed times are time marks of time_marks object passed at construction time with particular mask.
 *
 * The time step can be chosen either only for the next time level (this is default behavior of  next_time() ) or
 * for the largest possible time interval. The later possibility is necessary for explicit solvers, where the matrix used to perform one time step
 * depends on the time step and has to be modified, when the time step is changed.
 *
 * Information provided by time governor includes:
 * - actual time, last time
 * - actual time step
 * - number of the time level
 * - time comparison
 *
 * TODO: should be initialized directly from JSON input
 *
 */

class TimeMarks;

class TimeGovernor
{
public:
    /**
     * Constructor - constructor for unsteady solvers
     *
     * @param init_time - initial time (the solution stands on the initial value before this time)
     * @param end_time - final time of the particular equation
     * @param marks - reference to TimeMarks object which will be used for fixed times
     * @param fixed_time_mask - TimeMark mask used to select fixed times from all time marks
     *
     */
   TimeGovernor(const double init_time, const  double end_time, TimeMarks &marks, const TimeMark::Type fixed_time_mask = TimeMark::strict);

   /**
    * Default constructor - only for steady solvers
    *
    * We can have "zero step" steady problem (no computation, e.g. EquationNothing) and one step steady problem
    * (e.g. steady water flow). This constructor is only shortcut for those cases which set TimeMarks to NULL, end time
    * to infinity and time step to infinity.
    *
    * Optionally you can set initial time for "one step" steady problems.
    * First call of next_time() push the actual time to infinity.
    *
    * However, you have to use full constructor for the "steady problem" that have time variable input data.
    */
   TimeGovernor(double init_time = inf_time);

   /**
    * Permanent constrain for time step.
    */
   void set_permanent_constrain( double min_dt, double max_dt);

    /**
     * Set upper constrain for the choice of the next time step.
     *
     * When the time value is determined, the constrain is set to infinity. Until the next call of next_time you can
     * set constrains. The minimum of them is used for choice of the next time step.
     */
    void set_constrain(double dt_constrain);

    /**
     * Getter for time marks.
     */
    inline TimeMarks &marks() const
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
     * Add a time that has to be meet by the time governor.
     * The time is inserted into priority queue of fixed times in which the solution has to be computed.
     */
    //void set_fix_time(double fix_time);

    /**
     * Add more equidistant fix times.
     */
    //void set_fix_times(double first_fix_time, double fix_interval);

    /**
     *  Set the overhead in execution time that the solver has to do
     *  when the time step is changed. @fn dt_change_overhead is the overhead relative to
     *  the execution time of computation of one time step. Time governor tries to optimize
     *  the overall execution time according to this overhead estimate.
     *
     *  The value has to be positive. Negative or zero values turns off algorithm that keeps
     *  time step for longer periods.
     *
     *  Change in the time step can be tested by @fn is_changed_dt() and end of the time interval with fixed dt is provided
     *  by @fn end_of_fixed_dt()
     */

    void set_dt_change_overhead(double overhead = 0.0)
        { dt_change_overhead= overhead > 0.0 ? overhead : 0.0;}

    inline double end_of_fixed_dt() const
        {return end_of_fixed_dt_interval;}

    inline bool is_changed_dt() const
        {return dt_changed;}

    /**
     * Choice of the next time step according to the dt constrains, fixed times and overhead.
     */
    void next_time();

    /**
     * End of actual time interval; i.e. where the solution is computed.
     */
    inline double t() const
        {return time;}

    /**
     * Beginning of the actual time interval; i.e. the time of last computed solution.
     */
    inline double last_t() const
        {return last_time;}

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

    /// End time.
    inline double end_time() const
    { return end_time_; }

    inline bool is_end() const
        {return this->ge(end_time_); }

    /**
     * Performs comparison time > other_time, i.e. time is strictly greater then given parameter other_time
     * with precision relative to the time step. Precision is one percent of the time step.
     */
    inline bool gt(double other_time) const
        { return time >= other_time + comparison_fracture(); }

    /**
     * Performs comparison time >= other_time, with precision relative to the time step. See @fn gt
     */
    inline bool ge(double other_time) const
        { bool result = time + comparison_fracture() >= other_time;
          //DBGMSG("ge: %d %f %f %f\n", result, time, other_time, time_step);
          return result;
        }

    /**
     * Performs comparison time < other_time, with precision relative to the time step. See @fn gt
     */
    inline bool lt(double other_time) const
        { return time + comparison_fracture() <= other_time; }

    /**
     * Performs comparison time <= other_time, with precision relative to the time step. See @fn gt
     */
    inline bool le(double other_time) const
        { return time <= other_time + comparison_fracture(); }

    /**
     * Returns the time level.
     */
    inline int tlevel() const
        {return time_level;}


    void view() const
    {
        DBGMSG(" level: %d end time: %f time: %f step: %f\n",time_level, end_time_, time, time_step);
    }

    /// Infinity time used for steady case.
    static const double inf_time;

private:
    inline double comparison_fracture() const
    {
        if (0 < time_step && time_step <=numeric_limits<double>::max() ) return comparison_precision * time_step;
        else return numeric_limits<double>::epsilon();
    }


    /// we consider time difference is zero if it is less then comparison_precision * time_step
    static const double comparison_precision;
    /// technical bound for the time step given by finite precision
    static const double time_step_lower_bound;

    /// Number of time_next calls, i.e. total number of performed time steps.
    int time_level;
    /// End of actual time interval; i.e. where the solution is computed.
    double time;
    /// Beginning of the actual time interval; i.e. the time of last computed solution.
    double last_time;
    /// End of interval if fixed time step. (defers from @var time only if overhead is positive)
    double end_of_fixed_dt_interval;
    /// End time of the simulation.
    double end_time_;

    /// Length of actual time interval; i.e. the actual time step.
    double time_step;
    /// time step just before last_time
    double last_time_step;
    /// changed dt flag
    bool dt_changed;

    /// Overhead when we change dt. Relative to the computation of one time step.
    double dt_change_overhead;

    /// Upper constrain for the choice of the next time step. Relaxed after every dt choice.
    double time_step_constrain;
    /// Permanent upper limit for the time step.
    double max_time_step;
    /// Permanent lower limit for the time step.
    double min_time_step;

    /**
     * When the next time is chosen we need only the lowest fix time. Therefore we use
     * minimum priority queue of doubles based on the vector container.
     */
    TimeMarks * const time_marks;
    const TimeMark::Type fixed_time_mark_mask;

};

#endif /* TIME_HH_ */
