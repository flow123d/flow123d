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
 * $Id: sparse_graph.hh 822 2010-12-17 15:18:04Z jiri.jenicek $
 * $Revision: 822 $
 * $LastChangedBy: jiri.jenicek $
 * $LastChangedDate: 2010-12-17 16:18:04 +0100 (Fri, 17 Dec 2010) $
 *
 * @file
 * @brief Basic time management class.
 */

#ifndef TIME_HH_
#define TIME_HH_

#include <algorithm>
#include <queue>
#include "system.hh"
/**
 * @brief
 * This class should provide basic time management functionality for unsteady solvers.
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
    TimeGovernor(double time_init, double dt, double end_t);

    /**
     * Set constrain for the choice of the next time step.
     *
     * When the time value is determined, the constrain is set to infinity. Until the next call of next_time you can
     * set constrains. The minimum of them is used for choice of the next time step.
     */
    void constrain_dt(double dt_constrain);

    /**
     * Add an time into priority queue of fixed times in which the solution has to be computed.
     */
    void set_fix_time(double fix_time);

    /**
     * Choice of the next time step according to the suggested dt and fixed times.
     */
    void next_time();

    /**
     * Returns the time value.
     */
    inline double t() const
        {return time;}

    /**
     * Returns the time step value.
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

    double time;
    double last_time;
    double end_time;

    double time_step;
    double time_step_constrain;
    double min_time_step;
    double max_time_step;
    /**
     * When the next time is chosen we need only the lowest fix time. Therefore we use
     * minimum priority queue of doubles based on the vector container.
     */
    std::priority_queue<double, vector<double>, greater<double> > fix_times;

};

#endif /* TIME_HH_ */
