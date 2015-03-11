/*
 * TimerData.hh
 *
 *  Created on: 2. 3. 2015
 *      Author: jan-hybs
 */

#ifndef TIMERDATA_HH_
#define TIMERDATA_HH_

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
 */
class TimerData {
    public:

        /**
         * Method returns time in seconds in double precision
         */
        double to_time (void) const;

        /**
         * Debug-only method returning double value from 'toTime' method as string
         */
        string to_string (void);

        /**
         * Static method which generates TimerData object and sets its ticks value to current time
         */
        static TimerData get_time (void);

        /**
         * Static method for initializing timers (necessary for Windows timer only)
         */
        static void init (void);

        /**
         * Overloaded operator for addition. Allows adding TimerData objects using '+' sign
         */
        TimerData operator+ (const TimerData &right);

        /**
         * Overloaded operator for subtraction. Allows subtracting TimerData objects using '-' sign
         */
        TimerData operator- (const TimerData &right);

    private:
        /**
         * Constructor will populate object with current time
         */
        TimerData ();

        /**
         * Method which returns current time ticks. Depending on used timer value may represent
         * nanoseconds or internal ticks / CPU frequency
         */
        long long get_ticks (void) const;

        /**
         * Setter for ticks variable
         */
        void set_ticks (long long ticks);

        /**
         * Static method which gets current ticks and converts them to long long int
         */
        static long long get_current_ticks ();

        /**
         * Static variable representing initialization
         */
        static bool inited;
        /**
         * Internal variable for storing actual ticks
         */
        long long ticks_;

#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
        /**
         * Variable for storing CPU Frequency
         */
        static LARGE_INTEGER frequency;

        /**
         * Init function which sets current CPU frequency
         * and returns it
         */
        static LARGE_INTEGER get_frequency;
#endif //FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER

};
#endif /* TIMERDATA_HH_ */
