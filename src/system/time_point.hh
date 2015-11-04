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
 * @file    time_point.hh
 * @brief   
 */

#ifndef TIMEPOINT_HH_
#define TIMEPOINT_HH_

#include "config.h"
#include <string>
#include <time.h>

/**
 * Include either Windows or chrono lib
 */
#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
#include <windows.h>
#else
#include <chrono>
#endif //FLOW123D_HAVE_TIMER_CHRONO_HIGH_RESOLUTION

using namespace std;

/**
 * Class TimerData serves for getting current time in ticks units.
 * TimerData overloads standard + and - operators so object can be
 * subtracted/added from/to one another.
 *
 * TimerData uses one if the following timers based
 * on set definition namely:
 *
 * 1) FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
 * It is used Timer for Windows platform only.
 *
 * Method QueryPerformanceCounter populates given
 * LARGE_INTEGER object with timer informations.
 * Resolution if this timer is several hundreds ns.
 * Tests have shown average about 300 ns.
 * Maximum resolution is 1 ns.
 *
 *
 * 2) FLOW123D_HAVE_TIMER_CHRONO_HIGH_RESOLUTION (default without definition set)
 * It is used standard multiplatform Chrono Timer
 * which is available since C++11 (compile flag -std=cxx11 required!)
 *
 * Resolution of timer depends on HW but tests shown results fluctuating around 1000 ns.
 * Maximum resolution is 1 ns.
 *
 *
 * Note:
 * Standard clock() method available in c++ yields
 * resolution in units of milliseconds (roughly 10 ms)
 * but some compilers replaces during compilation
 * method clock() with better Timer having resolution 1000 ns.
 * Meaning clock() method is not sufficient or stable
 * for time measurements.
 * Maximum resolution is 1000 ns (limited by value of CLOCKS_PER_SEC).
 *
 * gettimeofday
 * The actual resolution of gettimeofday() depends on the hardware architecture. Intel
 * processors as well as SPARC machines offer high resolution timers that measure microseconds.
 * Other hardware architectures fall back to the system’s timer, which is typically set to 100 Hz.
 * In such cases, the time resolution will be less accurate. gettimeofday() can result in incorrect
 * timings if there are processes on your system that change the timer
 * Possible solution may be other clock
 *
 * clock_gettime (CLOCK_MONOTONIC) with clock_getres() method
 * but clock_gettime is present only on newest Linux. other system have only gettimeofday()
 */
class TimePoint {
    public:

        /**
         * Constructor will populate object with current time
         */
        TimePoint ();

        /**
         * Overloaded operator for subtraction. Allows subtracting TimerPoint objects using '-' sign
         * Used for determining interval between two TimerPoints
         *
         * Returns duration in seconds
         */
        double operator- (const TimePoint &right);

        /**
         * Internal variable for storing actual ticks
         */
        long long ticks;

    private:

#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
        /**
         * Variable for storing CPU Frequency
         */
        static LARGE_INTEGER frequency;

        /**
         * Init function which sets current CPU frequency
         * and returns it
         */
        static LARGE_INTEGER get_frequency ();
#endif //FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER

};
#endif /* TIMEPOINT_HH_ */
