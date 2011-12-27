/* 
 * File:   time..hh
 * Author: jb
 *
 * Created on May 29, 2010, 8:40 PM
 */

#ifndef _TIME__HH
#define	_TIME__HH

#include <queue>
#include "base/parameter_handler.h"


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
    SolverTime(ParameterHandler &prm)
                :step_number(0),
                suggested_dt(-1.0)


    {

        start_time = prm.get_double("t_init");
        time = start_time;
        end_time = prm.get_double("t_end");
        time_step = prm.get_double("dt_init");
        step_min = prm.get_double("dt_min");
        step_max = prm.get_double("dt_max");

        target_times.push(end_time);
    }

    void add_target_time(double target_time)
    {
        target_times.push(target_time);
    }

    bool inc()
    {
        bool is_end = time > end_time * (1 - 0.0000001);
        if (is_end) return false;

        if (suggested_dt > 0.0) time_step=suggested_dt;
        suggested_dt=-1.0;

        while (target_times.top() < time) target_times.pop();

        time_step=(target_times.top() - time)/ceil( (target_times.top() - time)/time_step );
        if (time_step < step_min) time_step = step_min;
        if (time_step > step_max) time_step = step_max;

        time+=time_step;
        step_number++;

        return true;
    }

    bool reinc_time(double factor)
    {
        bool result =true;
        time-=time_step;
        time_step=time_step*factor;
        if (time_step < step_min) {time_step = step_min; result=false;}
        if (time_step > step_max) time_step = step_max;

        time+=time_step;
        return result;
    }

    void scale_time_step(double factor) 
    {suggested_dt=time_step*factor;}


    double t() { return time;}
    double dt() {return time_step;}
    int n_step() {return step_number;}
    double end_t() {return end_time;}
};

#endif	/* _TIME__HH */

