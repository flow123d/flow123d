/*
 * TimerData.cc
 *
 *  Created on: 2. 3. 2015
 *      Author: jan-hybs
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>

#include "config.h"
#include "timer_data.hh"

using namespace std;


#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER

    #include <windows.h>

    TimePoint::TimePoint () {
        LARGE_INTEGER time;
        QueryPerformanceCounter (&time);
        this->ticks = time.QuadPart;
    }


    // initialize static frequency when using QueryPerformanceCounter
    LARGE_INTEGER TimePoint::get_frequency () {
        LARGE_INTEGER frequency;
        QueryPerformanceFrequency(&frequency);
        return frequency;
    }
    LARGE_INTEGER TimePoint::frequency = TimePoint::get_frequency ();


    double TimePoint::operator- (const TimePoint &right) {
        double difference = this->ticks - right.ticks;
        return difference / (TimePoint::frequency.QuadPart);
    }


#else

    #include <chrono>


    TimePoint::TimePoint () {
        chrono::time_point<chrono::high_resolution_clock> time = chrono::high_resolution_clock::now ();
        this->ticks = chrono::duration_cast<std::chrono::nanoseconds> (time.time_since_epoch ()).count ();
    }


    double TimePoint::operator- (const TimePoint &right) {
        double difference = this->ticks - right.ticks;
        return difference / (1 * 1000 * 1000 * 1000);
    }

#endif //FLOW123D_HAVE_TIMER_CHRONO_HIGH_RESOLUTION
