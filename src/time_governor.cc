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

#include <time_governor.hh>
#include <algorithm>

/*
 * TODO:
 * TimeGovernor should be constructed from JSON object.
 */
TimeGovernor::TimeGovernor(double time_init, double dt, double end_t)
{
    INPUT_CHECK( DBL_GT(dt, 0.0),"Time step has to be greater than ZERO\n");
    time=time_init;
    end_time=end_t;

    time_step=dt;
    min_time_step=dt;
    max_time_step=dt;
    time_step_constrain = min(end_time-time, max_time_step);

    time_level=0;
    fix_times.push(end_t);
}

void TimeGovernor::constrain_dt(double dt_constrain)
{
    time_step_constrain = min(time_step_constrain, dt_constrain);
}

void TimeGovernor::set_fix_time(double fix_time)
{
    fix_times.push(end_time);
}

void TimeGovernor::next_time()
{
    if (is_end()) return;
    last_time=time;

    // move to the next fix time
    while ( this->ge(fix_times.top()) ) fix_times.pop();

    // compute step to next fix time and apply constrains
    double full_step = fix_times.top() - last_time;
    time_step = min(full_step, time_step_constrain);
    time_step = min(time_step, max_time_step);
    time_step = max(time_step,min_time_step);

    // round the time step to have integer number of steps till next fix time
    // this always select shorter time step
    int n_steps = ceil( full_step / time_step );
    time_step = full_step / n_steps;
    time += time_step;

    // reset time step constrain
    time_step_constrain = min(end_time-time, max_time_step);
}
