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
    return ( stream << mark.time()<<": 0x" << hex << mark.mark_type().bitmap_ << dec );
}

const TimeMark::Type TimeMark::every_type =  TimeMark::Type(~0x0, std::numeric_limits<unsigned char>::min() );
const TimeMark::Type TimeMark::none_type =  TimeMark::Type(0x0, std::numeric_limits<unsigned char>::max() );


// ------------------------------------------------------
// implementation of members of class TimeMarks
// ------------------------------------------------------

TimeMarks::TimeMarks()
{
    this->reinit();
}



void TimeMarks::reinit()
{
	marks_.resize(1);
    marks_[0].clear();
    next_mark_type_ = TimeMark::Type(0x00, 1);

    // add predefined base mark types
    type_fixed_time_ = TimeMark::Type(0x01, 0);
    type_output_ = TimeMark::Type(0x02, 0);
    type_input_ = TimeMark::Type(0x04, 0);
    type_balance_output_ = TimeMark::Type(0x08, 0);

    // insert start and end stoppers
    marks_[0].push_back(TimeMark(-INFINITY, TimeMark::every_type));
    marks_[0].push_back(TimeMark(+INFINITY, TimeMark::every_type));
}

TimeMark::Type TimeMarks::new_mark_type() {
	ASSERT(next_mark_type_.equation_index_ != 0).error("Can not allocate new mark type. The limit is 255 mark types.\n");
    TimeMark::Type current_type = next_mark_type_;

    // add new subvector to marks_, insert start and end stoppers
    marks_.push_back( std::vector<TimeMark>() );
    marks_[next_mark_type_.equation_index_].push_back(TimeMark(-INFINITY, TimeMark::Type(~0x0, next_mark_type_.equation_index_ )));
    marks_[next_mark_type_.equation_index_].push_back(TimeMark(+INFINITY, TimeMark::Type(~0x0, next_mark_type_.equation_index_ )));

    next_mark_type_.equation_index_++;
    return current_type;
}

TimeMark TimeMarks::add(const TimeMark &mark) {
    // find first mark with time greater or equal to the new mark
    vector<TimeMark>::iterator first_ge = std::lower_bound(marks_[mark.mark_type_.equation_index_].begin(),
                                                           marks_[mark.mark_type_.equation_index_].end(), mark);

    // check equivalence with found mark
    if (fabs(first_ge->time() - mark.time()) < TimeGovernor::time_step_precision) {
	//if times are "equal" do bitwise OR with the mark type at the first_ge iterator position
        first_ge->add_to_type(mark.mark_type());
        return *first_ge;
    }
    // possibly check equivalence with previous mark
    if (first_ge != marks_[mark.mark_type_.equation_index_].begin()) {
        vector<TimeMark>::iterator previous = first_ge;
        --previous;
        if (fabs(previous->time() - mark.time()) < TimeGovernor::time_step_precision) {
            previous->add_to_type(mark.mark_type());
            return *previous;
        }
    }

    marks_[mark.mark_type_.equation_index_].insert(first_ge, mark);
    return mark;
}

void TimeMarks::add_time_marks(double time, double dt, double end_time, TimeMark::Type type) {
	ASSERT(end_time != TimeGovernor::inf_time).error("Can not add time marks on infinite interval.\n");
	ASSERT_GT(dt, numeric_limits<double>::epsilon()).error("TimeMark's step less then machine precision.\n");

	unsigned int n_steps=((end_time-time)/dt + TimeGovernor::time_step_precision);
	for (unsigned int i = 0; i<=n_steps;i++) {
		auto mark = TimeMark(time+i*dt, type);
		add(mark);
	}
}


void  TimeMarks::add_to_type_all(TimeMark::Type filter_type, TimeMark::Type add_type) {
    ASSERT(filter_type.equation_index_ != 0).error();
    ASSERT( (filter_type.equation_index_ == add_type.equation_index_) || (add_type.equation_index_ == 0) )
            ((unsigned int)filter_type.equation_index_)((unsigned int)add_type.equation_index_).error();

	for(auto &mark : marks_[filter_type.equation_index_])
        if (mark.match_mask(filter_type)) mark.add_to_type(add_type);

}

/*
bool TimeMarks::is_current(const TimeStep &time_step, const TimeMark::Type &mask) const
{
    return ( current(time_step, mask) != this->end(mask) );

    if (tg.t() == TimeGovernor::inf_time) return tg.is_end();
    const TimeMark &tm = *last(tg, mask);

    return tg.step().lt(tm.time() + tg.dt()); // last_t + dt < mark_t + dt

}*/


TimeMarks::iterator TimeMarks::current(const TimeStep &time_step, const TimeMark::Type &mask) const
{
    //if (time_step.end() == TimeGovernor::inf_time) return tg.is_end();
    auto it = last(time_step, mask);
    if ( time_step.contains(it->time()) ) return it;
    else return this->end(mask);
}

TimeMarks::iterator TimeMarks::next(const TimeGovernor &tg, const TimeMark::Type &mask) const
{
	ASSERT(tg.equation_mark_type().equation_index_ == mask.equation_index_)
			((unsigned int)tg.equation_mark_type().equation_index_)((unsigned int)mask.equation_index_).error();

	unsigned char eq_index = mask.equation_index_;
    // first time mark which does not compare less then then actual tg time
    vector<TimeMark>::const_iterator first_ge =
            std::lower_bound(marks_[eq_index].begin(), marks_[eq_index].end(), TimeMark(tg.t(),mask));
    while (  ! tg.step().lt(first_ge->time()) || ! first_ge->match_mask(mask) ) {
        ++first_ge;
        if (first_ge == marks_[eq_index].end()) { --first_ge; break; }
    }
    return TimeMarksIterator(marks_[eq_index], first_ge, mask);
}

TimeMarks::iterator TimeMarks::last(const TimeStep &time_step, const TimeMark::Type &mask) const
{
	unsigned char eq_index = mask.equation_index_;
    // first time mark which does compare strictly greater then actual tg time
    vector<TimeMark>::const_iterator first_ge =
            std::lower_bound(marks_[eq_index].begin(), marks_[eq_index].end(), TimeMark(time_step.end()+0.01*time_step.length(),mask));
    while ( ! time_step.ge(first_ge->time()) || ! first_ge->match_mask(mask) ) {
        --first_ge;
    }
    return TimeMarksIterator(marks_[eq_index], first_ge, mask);
}

TimeMarks::iterator TimeMarks::last(const TimeGovernor &tg, const TimeMark::Type &mask) const
{
    return last(tg.step(), mask);
}


TimeMarks::iterator TimeMarks::last(const TimeMark::Type &mask) const
{
	auto it = TimeMarksIterator(marks_[mask.equation_index_], --marks_[mask.equation_index_].end(), mask); // +INF time mark
	--it;
	return it;
}



TimeMarks::iterator TimeMarks::begin(TimeMark::Type mask) const
{
	return TimeMarksIterator(marks_[mask.equation_index_], marks_[mask.equation_index_].begin(), mask);
}



TimeMarks::iterator TimeMarks::begin(unsigned char ei) const
{
	return this->begin( TimeMark::Type(TimeMark::every_type.bitmap_, ei) );
}



TimeMarks::iterator TimeMarks::end(TimeMark::Type mask) const
{
	return TimeMarksIterator(marks_[mask.equation_index_], --marks_[mask.equation_index_].end(), mask);
}



TimeMarks::iterator TimeMarks::end(unsigned char ei) const
{
	return this->end( TimeMark::Type(TimeMark::every_type.bitmap_, ei) );
}



ostream& operator<<(ostream& stream, const TimeMarks &marks)
{
	stream << "==============================" << endl;
    stream << "time marks:" << endl;
    for(unsigned int i=0; i<marks.marks_.size(); ++i) {
    	stream << "------------------------------" << endl;
    	stream << "equation_index: " << i << endl;
    	for(vector<TimeMark>::const_iterator it = marks.marks_[i].begin(); it != marks.marks_[i].end(); ++it)
			stream << *it << endl;
    }
	stream << "==============================" << endl;
    return stream;
}
