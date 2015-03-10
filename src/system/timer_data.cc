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

#include "config.h"
#include "timer_data.hh"

using namespace std;

#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
#include <windows.h>
#else
#include <chrono>
#endif //FLOW123D_HAVE_TIMER_CHRONO_HIGH_RESOLUTION

bool TimerData::inited = false;

/*
 * Constructor
 */

TimerData::TimerData () :
        ticks_ (TimerData::get_current_ticks ()) {
}

/*
 * Public methods
 */

double TimerData::to_time (void) const {
#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
    return ((double)this->get_ticks()) / (TimerData::frequency.QuadPart);
#else
    return ((double) this->get_ticks ()) / (1 * 1000 * 1000 * 1000);
#endif //FLOW123D_HAVE_TIMER_CHRONO_HIGH_RESOLUTION
}

string TimerData::to_string (void) {
    char buffer[50];
    sprintf (buffer, "%1.9f", this->to_time ());
    return buffer;
}

TimerData TimerData::operator+ (const TimerData &right) {
    TimerData result;
    result.set_ticks (this->ticks_ + right.ticks_);
    return result;
}

TimerData TimerData::operator- (const TimerData &right) {
    TimerData result;
    result.set_ticks (this->ticks_ - right.ticks_);
    return result;
}

/*
 * Static methods
 */

TimerData TimerData::get_time () {
    TimerData result;
    // ticks are already set in constructor
    return result;
}
void TimerData::init () {
    if (TimerData::inited) return;
    // init routine
    TimerData::inited = true;
}

/*
 * Private methods
 */

long long TimerData::get_ticks () const {
    return this->ticks_;
}

void TimerData::set_ticks (long long ticks) {
    this->ticks_ = ticks;
}

long long TimerData::get_current_ticks () {
#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
    LARGE_INTEGER time;
    QueryPerformanceCounter (&time);
    return time.QuadPart;
#else
    chrono::time_point<chrono::high_resolution_clock> time = chrono::high_resolution_clock::now ();
    return chrono::duration_cast<std::chrono::nanoseconds> (time.time_since_epoch ()).count ();
#endif //FLOW123D_HAVE_TIMER_CHRONO_HIGH_RESOLUTION
}

#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
LARGE_INTEGER TimerData::get_frequency () {
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    return frequency;
}
// initialize frequency field
LARGE_INTEGER TimerData::frequency = TimerData::get_frequency ();
#endif //FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER

std::ostream& operator<< (std::ostream &strm, TimerData &right) {
    return strm << right.to_string ();
}

