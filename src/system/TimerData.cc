/*
 * TimerData.cc
 *
 *  Created on: 2. 3. 2015
 *      Author: jan-hybs
 */
#include "TimerData.hh"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using namespace std;

#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
#include <windows.h>
LARGE_INTEGER TimerData::Frequency;
#else
#include <chrono>
#endif //FLOW123D_HAVE_TIMER_CHRONO_HIGH_RESOLUTION

bool TimerData::inited = false;

/*
 * Constructor
 */

TimerData::TimerData () :
        ticks_ (TimerData::getCurrentTicks ()) {
}

/*
 * Public methods
 */

double TimerData::toTime (void) const {
#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
    return ((double)this->getTicks()) / (TimerData::Frequency.QuadPart);
#else
    return ((double) this->getTicks ()) / (1 * 1000 * 1000 * 1000);
#endif //FLOW123D_HAVE_TIMER_CHRONO_HIGH_RESOLUTION
}

string TimerData::toString (void) {
    char buffer[50];
    sprintf (buffer, "%1.9f", this->toTime ());
    return buffer;
}

TimerData TimerData::operator+ (const TimerData &right) {
    TimerData result;
    result.setTicks (this->ticks_ + right.ticks_);
    return result;
}

TimerData TimerData::operator- (const TimerData &right) {
    TimerData result;
    result.setTicks (this->ticks_ - right.ticks_);
    return result;
}

/*
 * Static methods
 */

TimerData TimerData::getTime () {
    TimerData result;
    // ticks are already set in constructor
    return result;
}
void TimerData::init () {
    if (TimerData::inited) return;

    // set Frequency on Windows platform
#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
    QueryPerformanceFrequency(&TimerData::Frequency);
//        cout << "Windows used" << endl;
#else
//        cout << "Chrono used" << endl;
#endif //FLOW123D_HAVE_TIMER_CHRONO_HIGH_RESOLUTION
    TimerData::inited = true;
}

/*
 * Private methods
 */

long long TimerData::getTicks () const {
    return this->ticks_;
}

void TimerData::setTicks (long long ticks) {
    this->ticks_ = ticks;
}

long long TimerData::getCurrentTicks () {
#ifdef FLOW123D_HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
    LARGE_INTEGER time;
    QueryPerformanceCounter (&time);
    return time.QuadPart;
#else
    chrono::time_point<chrono::high_resolution_clock> time = chrono::high_resolution_clock::now ();
    return chrono::duration_cast<std::chrono::nanoseconds> (time.time_since_epoch ()).count ();
#endif //FLOW123D_HAVE_TIMER_CHRONO_HIGH_RESOLUTION
}

std::ostream& operator<< (std::ostream &strm, TimerData &right) {
    return strm << 564;
}

