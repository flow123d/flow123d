/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    time_point.cc
 * @brief   
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>

#include "config.h"
#include "time_point.hh"

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
