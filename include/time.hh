/* 
 * File:   time..hh
 * Author: jb
 *
 * Created on May 29, 2010, 8:40 PM
 */

#ifndef _TIME__HH
#define	_TIME__HH

#include <queue>


class SolverTime {
    // parameters
    double start_time;
    double end_time;

    // variables
    double time;
    double time_step;
    double step_min, step_max;
    unsigned int step_number;
    std::priority_queue<double> target_times;
    double suggested_dt;

public:
    SolverTime(double start, double end, double dt_init)
            :   start_time(start),
                end_time(end),
                time(start),
                time_step(dt_init),
                step_number(0),
                suggested_dt(-1.0)

    { target_times.push(end); }

    void add_target_time(double target_time)
    {
        target_times.push(target_time);
    }

    bool inc()
    {
        if (suggested_dt > 0.0) time_step=suggested_dt;
        suggested_dt=-1.0;

        while (target_times.top() < time) target_times.pop();

        time_step=(target_times.top() - time)/ceil( (target_times.top() - time)/time_step );
        time+=time_step;
        step_number++;

        return (time <= end_time);
    }

    void reinc_time(double factor)
    {
        time-=time_step;
        time_step=time_step*factor;
        time+=time_step;
    }

    void scale_time_step(double factor) 
    {suggested_dt=time_step*factor;}


    double t() { return time;}
    double dt() {return time_step;}
    int n_step() {return step_number;}
};

#endif	/* _TIME__HH */

