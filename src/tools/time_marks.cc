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
 * @file    time_marks.cc
 * @brief   
 * @author  Jan Brezina
 */

#include <algorithm>
#include <limits>
#include "system/system.hh"
#include "system/global_defs.h"
#include "time_governor.hh"
#include "time_marks.hh"

// ------------------------------------------------------
// implementation of members of class TimeMark
// ------------------------------------------------------

ostream& operator<<(ostream& stream, const TimeMark &mark)
{
    //return ( stream << mark.time()<<": 0o" << oct << mark.mark_type() << dec ); //octal output
    return ( stream << mark.time()<<": 0x" << hex << mark.mark_type() << dec );
}

const TimeMark::Type TimeMark::every_type =  ~0x0;
const TimeMark::Type TimeMark::none_type =  0x0;


// ------------------------------------------------------
// implementation of members of class TimeMarks
// ------------------------------------------------------

TimeMarks::TimeMarks()
{
    this->reinit();
}



void TimeMarks::reinit()
{
    marks_.clear();
    next_mark_type_ = 0x1;

    // add predefined base mark types
    type_fixed_time_ = new_mark_type();
    type_output_ = new_mark_type();
    type_input_ = new_mark_type();

    // insert start and end stoppers
    marks_.push_back(TimeMark(-INFINITY, TimeMark::every_type));
    marks_.push_back(TimeMark(+INFINITY, TimeMark::every_type));
}

TimeMark::Type TimeMarks::new_mark_type() {
    ASSERT(next_mark_type_ != 0, "Can not allocate new mark type. The limit is 32 mark types.\n");
    TimeMark::Type current_type = next_mark_type_;

    next_mark_type_ <<= 1;
    return current_type;
}

void TimeMarks::add(const TimeMark &mark) {
    // find first mark with time greater or equal to the new mark
    vector<TimeMark>::iterator first_ge = std::lower_bound(marks_.begin(), marks_.end(), mark);

    // check equivalence with found mark
    if (fabs(first_ge->time() - mark.time()) < TimeGovernor::time_step_precision) {
	//if "equal" does bitwise OR with the mark type at the first_ge iterator position
        first_ge->add_to_type(mark.mark_type());
        return;
    }
    // possibly check equivalence with previous mark
    if (first_ge != marks_.begin()) {
        vector<TimeMark>::iterator previous = first_ge;
        --previous;
        if (fabs(previous->time() - mark.time()) < TimeGovernor::time_step_precision) {
            previous->add_to_type(mark.mark_type());
            return;
        }
    }

    marks_.insert(first_ge, mark);
}

void TimeMarks::add_time_marks(double time, double dt, double end_time, TimeMark::Type type) {
	ASSERT(end_time != TimeGovernor::inf_time, "Can not add time marks on infinite interval.\n");
	ASSERT(dt > numeric_limits<double>::epsilon(), "TimeMark's step less then machine precision.\n");

	unsigned int n_steps=((end_time-time)/dt + TimeGovernor::time_step_precision);
	for (unsigned int i = 0; i<=n_steps;i++) {
		auto mark = TimeMark(time+i*dt, type);
		add(mark);
	}
}

bool TimeMarks::is_current(const TimeGovernor &tg, const TimeMark::Type &mask) const
{
    if (tg.t() == TimeGovernor::inf_time) return tg.is_end();
    const TimeMark &tm = *last(tg, mask);

    return tg.step().lt(tm.time() + tg.dt()); // last_t + dt < mark_t + dt
}

TimeMarks::iterator TimeMarks::next(const TimeGovernor &tg, const TimeMark::Type &mask) const
{
    // first time mark which does not compare less then then actual tg time
    vector<TimeMark>::const_iterator first_ge =
            std::lower_bound(marks_.begin(), marks_.end(), TimeMark(tg.t(),mask));
    while (  ! tg.step().lt(first_ge->time()) || ! first_ge->match_mask(mask) ) {
        ++first_ge;
    }
    return TimeMarksIterator(marks_, first_ge, mask);
}

TimeMarks::iterator TimeMarks::last(const TimeGovernor &tg, const TimeMark::Type &mask) const
{
    // first time mark which does compare strictly greater then actual tg time
    vector<TimeMark>::const_iterator first_ge =
            std::lower_bound(marks_.begin(), marks_.end(), TimeMark(tg.t()+0.01*tg.dt(),mask));
    while ( ! tg.step().ge(first_ge->time()) || ! first_ge->match_mask(mask) ) {
        --first_ge;
    }
    return TimeMarksIterator(marks_, first_ge, mask);
}



TimeMarks::iterator TimeMarks::last(const TimeMark::Type &mask) const
{
	auto it = TimeMarksIterator(marks_, --marks_.end(), mask); // +INF time mark
	--it;
	return it;
}



TimeMarks::iterator TimeMarks::begin(TimeMark::Type mask) const
{
	return TimeMarksIterator(marks_, marks_.begin(), mask);
}



TimeMarks::iterator TimeMarks::end(TimeMark::Type mask) const
{
	return TimeMarksIterator(marks_, --marks_.end(), mask);
}



ostream& operator<<(ostream& stream, const TimeMarks &marks)
{
    stream << "time marks:" << endl;
    for(vector<TimeMark>::const_iterator it = marks.marks_.begin(); it != marks.marks_.end(); ++it)
        stream << *it << endl;
    return stream;
}
