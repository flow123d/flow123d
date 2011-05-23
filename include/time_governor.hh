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

#include <algorithm>
#include <queue>
#include "system/system.hh"
/**
 * @brief
 * Basic time management functionality for unsteady solvers.
 *
 * This class provides only algorithms for selecting next time step, and information about current time step frame.
 * In particular fixed times are not time events, but can be added in order to meet the time of a time event. Time events are
 * managed by class TimeEvents.
 *
 * The time step can be chosen either only for the next time level (this is default behavior od @fn next_time) or
 * for the largest possible time interval. The later possibility is necessary for explicit solvers, where the matrix used to perform one time step
 * depends on the time step and has to be modified, when the time step is changed.
 *
 * This includes:
 * - actual time
 * - actual time step
 * - choice of next time step
 * - fixed times that has to be meet
 * - time level
 * - time comparison
 * TODO: should be initialized directly from JSON input
 *
 * Ideas: Distinguish init_time and start_time. The former one is initial value of the time value, the latter one is the beginning of the time frame
 * where an unsteady solver takes a place. Could be useful for short time processes
 *
 */

class TimeGovernor
{
public:
    /**
     * Constructor - from given values, set fixed time step.
     */
   TimeGovernor(double time_init, double min_dt, double max_dt, double end_t);

    /**
     * Set upper constrain for the choice of the next time step.
     *
     * When the time value is determined, the constrain is set to infinity. Until the next call of next_time you can
     * set constrains. The minimum of them is used for choice of the next time step.
     */
    void constrain_dt(double dt_constrain);

    /**
     * Add a time that has to be meet by the time governor.
     * The time is inserted into priority queue of fixed times in which the solution has to be computed.
     */
    void set_fix_time(double fix_time);

    /**
     * Add more equidistant fix times.
     */
    void set_fix_times(double first_fix_time, double fix_interval);

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
     * Length of actual time interval; i.e. the actual time step.
     */
    inline double dt() const
        {return time_step;}

    inline bool is_end() const
        {return this->ge(end_time); }

    /**
     * Performs comparison time > other_time, i.e. time is strictly greater then given parameter other_time
     * with precision relative to the time step. Precision is one percent of the time step.
     */
    inline bool gt(double other_time) const
        { return time > other_time + comparison_precision * time_step; }

    /**
     * Performs comparison time >= other_time, with precision relative to the time step. See @fn gt
     */
    inline bool ge(double other_time) const
        { return time > other_time - comparison_precision * time_step; }

    /**
     * Performs comparison time < other_time, with precision relative to the time step. See @fn gt
     */
    inline bool lt(double other_time) const
        { return time < other_time - comparison_precision * time_step; }

    /**
     * Performs comparison time <= other_time, with precision relative to the time step. See @fn gt
     */
    inline bool le(double other_time) const
        { return time < other_time + comparison_precision * time_step; }

    /**
     * Returns the time level.
     */
    inline int tlevel() const
        {return time_level;}


    void view() const
    {
        DBGMSG(" level: %d end time: %f time: %f step: %f\n",time_level, end_time, time, time_step);
    }

private:
    static const double comparison_precision=0.01;
    int time_level;

    /// End of actual time interval; i.e. where the solution is computed.
    double time;
    /// Beginning of the actual time interval; i.e. the time of last computed solution.
    double last_time;
    /// End of interval if fixed time step. (defers from @var time only if overhead is positive)
    double end_of_fixed_dt_interval;
    /// End time of the simulation.
    double end_time;

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
    std::priority_queue<double, vector<double>, greater<double> > fix_times;

};

#endif /* TIME_HH_ */
