/*
 * time.hh
 *
 *  Created on: Jan 19, 2011
 *      Author: jb
 */

#ifndef TIME_HH_
#define TIME_HH_

#include "system.hh"
/**
 * This class SHOULD provide basic time management functionality for unsteady solvers.
 *
 * This includes:
 * - actual time
 * - actual time step
 * - choice of next time step
 * - fixed times that has to be meet
 * - time level
 * - time comparison
 * - should be initialized directly from JSON input
 */

class TimeGovernor
{
public:
    TimeGovernor() : time(0.0), time_step(0.0), end_time(0.0), time_level(0.0) {}

    void setup(double time_init, double dt, double end_t)
        {
        INPUT_CHECK( DBL_GT(dt, 0.0),"Time step has to be greater than ZERO\n");
        time=time_init;
        time_step=dt;
        end_time=end_t;
        time_level=0;
        }

    inline double t()
        {return time;}

    inline double dt()
        {return time_step;}

    inline bool is_end()
        {return is_time(end_time); }

    /// Returns true if given time parameter is greater then
    /// actual time of governer minus one percent of the time step
    inline bool is_time(double some_time)
        {
            return ( (time-some_time) > -time_step*0.01);
        }

    int tlevel()
        {return time_level;}

    void next_time()
        {
            if (! is_end()) {
                time+=time_step; time_level++;
            }
        }

    void view()
    {
        DBGMSG(" level: %d end time: %f time: %f step: %f\n",time_level, end_time, time, time_step);
    }
private:
    int time_level;
    double time;
    double time_step;
    double end_time;
};

#endif /* TIME_HH_ */
