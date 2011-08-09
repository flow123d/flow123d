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
 * @ingroup application
 * @brief Basic time management class.
 */

#include "system/system.hh"
#include <time_governor.hh>
#include <time_marks.hh>

#include <limits>


const double TimeGovernor::comparison_precision = 0.01;
const double TimeGovernor::time_step_lower_bound = DBL_EPSILON;
const double TimeGovernor::inf_time =  numeric_limits<double>::infinity();

/*
 * TODO:
 * TimeGovernor should be constructed from JSON object.
 */
TimeGovernor::TimeGovernor(const double init_time,const  double end_time, TimeMarks &marks,const TimeMark::Type fixed_time_mask)
: last_time(init_time),
  end_time_(end_time),
  time_marks(&marks),
  fixed_time_mark_mask(fixed_time_mask)
{
    time=last_time;
    end_of_fixed_dt_interval=time;

    dt_changed=true;
    dt_change_overhead=-1.0; // turn off

    min_time_step=0;
    if (end_time_ != inf_time)  max_time_step=end_time_ - init_time;
    else max_time_step = inf_time;

    time_step=max_time_step;

    time_step_constrain = min(end_time_-time, max_time_step);

    time_level=0;
    time_marks->add( TimeMark(init_time, fixed_time_mark_mask) );
    time_marks->add( TimeMark(end_time_, fixed_time_mark_mask) );
}

TimeGovernor::TimeGovernor(double init_time)
: last_time(init_time),
  end_time_(inf_time),
  time(inf_time),
  time_marks(NULL),
  fixed_time_mark_mask(TimeMark::strict)

{
    end_of_fixed_dt_interval=time;

    dt_changed=true;
    dt_change_overhead=-1.0; // turn off

    min_time_step=0;
    if (end_time_ != inf_time)  max_time_step=end_time_ - init_time;
    else max_time_step = inf_time;

    time_step=max_time_step;
    time_step_constrain = min(end_time_-time, max_time_step);

    time_level=0;
}

void TimeGovernor::set_permanent_constrain( double min_dt, double max_dt)
{
    ASSERT( min_dt >= 0.0,"Minimal time step has to be greater than ZERO\n");
    ASSERT( max_dt >= min_dt,"Maximal time step has to be greater or equal to the minimal.\n");

    min_time_step=max(min_dt, time_step_lower_bound);
    max_time_step=min(max_dt, end_time_-time);
}

void TimeGovernor::set_constrain(double dt_constrain)
{
    time_step_constrain = min(time_step_constrain, dt_constrain);
}



/*
void TimeGovernor::set_fix_time(double fix_time)
{
    if (fix_time >= end_of_fixed_dt_interval)
        fix_times.push(fix_time);
    else
        xprintf(Warn, "Inserted fixed time %f less then end of fixed dt interval %f.\n",fix_time, end_of_fixed_dt_interval);
}
*/

/*
void TimeGovernor::set_fix_times(double first_fix_time, double fix_interval)
{
    for(double tt=first_fix_time; tt< end_time_; tt+=fix_interval) set_fix_time(tt);
}
*/

void TimeGovernor::next_time()
{
    DBGMSG("%f %f\n",time, end_time_);
    if (time == inf_time || is_end()) return;
    if (end_time_ == inf_time) {
        time = end_time_;
        return;
    }

    last_time=time;
    last_time_step = time_step;

    // jump to the first future fix time
    TimeMarks::iterator fix_time_it = time_marks->next(*this, fixed_time_mark_mask);

    // select algorithm for determination of time step
    if (dt_change_overhead <= 0.0) {
        // LOCAL DT CHOICE

        // compute step to next fix time and apply constrains
        double full_step = fix_time_it->time() - last_time;
        time_step = min(full_step, time_step_constrain);
        time_step = min(time_step, max_time_step);
        time_step = max(time_step,min_time_step);

        // round the time step to have integer number of steps till next fix time
        // this always select shorter time step
        int n_steps = ceil( full_step / time_step );
        time_step = full_step / n_steps;

        end_of_fixed_dt_interval=time;
        dt_changed= (last_time_step == time_step);
        time_step_constrain = min(end_time_-time, max_time_step); // reset time step constrain
    } else {
        // OVERHEAD OPTIMIZATION


        if (this->ge(end_of_fixed_dt_interval)) {
            // end of fixed time_step

            // simple solution that do not take overhead into account, only add fix points until they match the pattern
            //
            // possible problem: suppose fix_times: 0.5   1.0   2.0   3.0 ...
            // for reasonable overhead the optimal choice is 0.5 till time 1.0 and then time step 1.0
            // how to detect this? We should determine time_step for interval after proposed end_of_fixed_dt_interval
            // and compare total price per time for both intervals.
            //

              // compute step to next fix time and apply constrains
              double full_step = fix_time_it->time() - last_time;
              time_step = min(full_step, time_step_constrain);
              time_step = min(time_step, max_time_step);
              time_step = max(time_step,min_time_step);

              // round the time step to have integer number of steps till next fix time
              // this always select shorter time step
              int n_steps = ceil( full_step / time_step );
              time_step = full_step / n_steps;

              while ( fix_time_it != time_marks->end() &&
                      fabs( round(fix_time_it->time() / time_step) - fix_time_it->time()/ time_step ) <= comparison_precision ) ++fix_time_it;

              end_of_fixed_dt_interval=fix_time_it->time();
              dt_changed= (last_time_step == time_step);
              time_step_constrain = min(end_time_-time, max_time_step);         // reset time step constrain

            /*
             * following piece of code implements variant of Euclidead algorithm for GCD, but it appears that
             * this can not be used for our problem. Mathematical problem behind is:
             *
             * Find dt such that for every x_i there exists n_i such that | n_i * dt - x_i| < eps * dt
             * where x_i is given set of real numbers.
             *
             * Problem is that if the condition holds for some dt, it does not hold for dt/2 for example
             * 1) we have to check the condition for every proposed dt value
             * 2) the Eucleidean sequance probably do not produce god estimates of dt
             *
             */
/*
            // till the first fixed time we only shorten the time_step_constrain
            time_step=min(time_step_constrain, max_time_step);
            fixed_dt_interval = fix_times.top()-time;
            n_steps = ceil( fixed_dt_interval / time_step );
            time_step = fixed_dt_interval / n_steps;

            best_price = (dt_change_overhead + n_steps)/ fixed_dt_interval;
            end_of_fixed_dt_interval = fix_times.top();
            fixed_times.pop();

            // try further fixed points by greatest common interval divisor
            while (1) {
                fixed_dt_interval = fix_times.top()-time;
                try_dt=common_dt(fixed_dt_interval, time_step);
                n_steps = ceil( fixed_dt_interval / try_dt );
                price = (dt_change_overhead + n_steps)/ fixed_dt_interval;

                if (try_dt < min_time_step || price > best_price) break; // can not find good GCD

                // keep trying
                best_price=price;
                time_step = fixed_dt_interval / n_steps;
                end_of_fixed_dt_interval = fix_times.top();
                fix_times.pop();
            }

            // reset time step constrain
            time_step_constrain = min(end_time_-time, max_time_step); */
        } else {
            dt_changed= false;
        }


    }

    time+=time_step;
    time_level++;

}

