/*
 * time_governor_test.cpp
 *
 *  Created on: May 20, 2011
 *      Author: jb
 */

#include "system/system.hh"
#include "time_marks.hh"
#include "time_governor.hh"

#define EQUAL(a,b) INPUT_CHECK( (a) == (b), #a": %f and "#b":%f differs\n",a,b);

typedef int i3[3];

int main(int argc, char * argv[]) {
// test time marks
    TimeMarks tm;
    tm.add(TimeMark(-1.0,TimeMark::strict));
    tm.add(TimeMark(0.0,TimeMark::every_type));
    TimeMark::Type my_mark_type = tm.new_mark_type();
    TimeMark::Type your_mark_type = tm.new_strict_mark_type();
    tm.add_time_marks(1.0, 0.25, 2.0, my_mark_type);
    tm.add(TimeMark(3.0, my_mark_type  | TimeMark::strict));
    tm.add(TimeMark(4.0, your_mark_type | TimeMark::strict));

    cout << tm;

    TimeGovernor  tm_tg(&tm, -2.0, 10.0);

    TimeMarks::iterator it = tm.next(tm_tg, TimeMark::strict);
    EQUAL((it)->time(), -1.0);
    EQUAL((++it)->time(), 0.0);
    EQUAL((++it)->time(), 3.0);
    EQUAL((++it)->time(), 4.0);

    it = tm.next(tm_tg, your_mark_type);
    EQUAL(it->time(), 0.0);
    EQUAL((++it)->time(), 4.0);

    tm_tg.next_time();
    EQUAL(tm_tg.t(), -1.0);
    EQUAL(tm.is_current(tm_tg, my_mark_type), false);

    tm_tg.next_time();
    EQUAL(tm_tg.t(), 0.0);
    EQUAL(tm.is_current(tm_tg, my_mark_type), true);

    tm_tg.next_time();
    EQUAL(tm_tg.t(), 3.0);
    EQUAL(tm.is_current(tm_tg, my_mark_type), true);

    it = tm.last(tm_tg, my_mark_type);
    EQUAL(it->time(), 3.0);
    EQUAL((--it)->time(), 2.0);
    EQUAL((--it)->time(), 1.75);

// test time governor
    TimeGovernor  tg(&tm, 0.0, 10.0);
    tg.set_permanent_constrain(0.01,1.0);
    
    EQUAL(tg.t(), 0.0);
    tg.next_time();
    EQUAL(tg.dt(), 1.0);
    EQUAL(tg.t(), 1.0);
/*
    tg.set_fix_times(0.0, 0.5);
    tg.set_fix_time(1.9);

    tg.next_time();
    EQUAL(tg.t(), 1.5);
    tg.constrain_dt(0.3);
    tg.next_time();
    EQUAL(tg.t(),1.7);
    tg.constrain_dt(0.3);
    tg.next_time();
    EQUAL(tg.t(),1.9);
    tg.constrain_dt(0.3);
    tg.next_time();
    EQUAL(tg.t(),2.0);
    tg.constrain_dt(0.3);
    tg.next_time();
    EQUAL(tg.t(),2.25);

    tg.set_dt_change_overhead(100.0);
    tg.constrain_dt(0.3);
    tg.next_time();
    EQUAL(tg.dt(), 0.25);
    INPUT_CHECK(tg.is_changed_dt(), "dt do not change\n");
    EQUAL(tg.end_of_fixed_dt(), 10.0);
    EQUAL(tg.t(), 2.5);

    tg.constrain_dt(0.1);
    tg.next_time();
    EQUAL(tg.dt(), 0.25);
    EQUAL(tg.end_of_fixed_dt(), 10.0);
    INPUT_CHECK( !tg.is_changed_dt(), "dt does change\n");
    EQUAL(tg.t(), 2.75);

*/

}
