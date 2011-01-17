/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id: profiler.cc 842 2011-01-08 17:58:15Z tomas.bambuch $
 * $Revision: 842 $
 * $LastChangedBy: tomas.bambuch $
 * $LastChangedDate: 2011-01-08 18:58:15 +0100 (So, 08 led 2011) $
 *
 * @file
 * @brief  Profiler
 *
 */


#include <sys/param.h>

#include "profiler.hh"
#include "system.hh"


Profiler* Profiler::_instance = NULL;
map<string, TimerFrame*> TimerFrame::_frames;

Profiler::Profiler() {
    F_ENTRY;
    xprintf(Msg, "Profiler created  \n");

    actual_node = root = new Timer("", NULL);
    root->start(0);

    start_clock = clock();
}

Profiler::~Profiler() {
    if (root) {
        root->forced_end(get_time()); //stop all running timers
        xprintf(Msg, "-------------Profiler information-------------\n");
        root->print(0);
        xprintf(Msg, "----------------------------------------------\n");
    }
}

void Profiler::start(string tag) {
    F_ENTRY;

    /**
     * starts particular timing period:
     * - if the tag is new:
     *        make new instance of timing class and connect it to actual leaf of the profiling tree
     *        register into map
     *        if it is the first timing, set it as a root
     * - if tag exists, change actual leaf, increment call count of the timing object
     */

    ASSERT(actual_node, "No active timing node!");

    map<string, Timer*>::const_iterator i = tag_map.find(tag);

    if (i == tag_map.end()) { //the key was not found

        Timer *tmr = new Timer(tag, actual_node);
        tmr->start(get_time());
        tag_map.insert(std::make_pair(tag, tmr));

        actual_node->insert_child(tmr);

        actual_node = tmr;

    }
    else {
        i->second->start(get_time());

        actual_node = i->second;
    }
}

void Profiler::end(string tag) {

    /**
     * - check if tag match tag of actual timing object, if not print warning and proceed to matchin parent, while closing all open timings
     * - else close actual timing and proceed to the parent node in profiling tree
     */

    ASSERT(actual_node, "No active timing node!");

    if (tag == "") {
        //don't allow to close the root timer
        if (actual_node != root) {
            //close actual timing and all of its children
            if (actual_node->end(get_time())) {

                actual_node = actual_node->parent();
            }
        }
    }
    else if (actual_node->tag() == tag) {
        //close actual timing and all of its children
        if (actual_node->end(get_time())) {

            actual_node = actual_node->parent();
        }
    }
    else {
        map<string, Timer*>::const_iterator i = tag_map.find(tag);
        if (i != tag_map.end()) {
            if (i->second->end(get_time())) {
                actual_node = i->second->parent();
            }
        }
        else {
            //ERROR - the key has not been found
        }
    }

}

double Profiler::get_time() {
    return 1000 * ((double) (clock() - start_clock)) / CLOCKS_PER_SEC;
}

Timer::Timer(string tag, Timer* parent) {

    timer_tag = tag;
    running = false;
    parent_timer = parent;
    start_count = 0;
    start_time = 0;
    cumul_time = 0;
    count = 0;
}

void Timer::start(double time) {
    F_ENTRY;

    count++;
    start_count++;
    if (!running) {

        xprintf(Msg, "Starting timer %s \n", this->tag().c_str());

        start_time = time;
        running = true;
    }
}

bool Timer::end(double time) {
    //close the timer only if start() has been called as many times as end()

    start_count = MAX(start_count - 1, 0);
    if (start_count == 0) {
        if (running)
            xprintf(Msg, "Stopping timer %s \n", this->tag().c_str());
        this->stop(time);
        return true;
    }

    return false;
}

void Timer::forced_end(double time) {
    if (running)
        xprintf(Msg, "Forced stopping timer %s \n", this->tag().c_str());
    start_count = 0;
    this->stop(time);
}

void Timer::stop(double time) {
    if (running && start_count == 0) {
        running = false;
        cumul_time += time - start_time;
        //force to close all of its children
        //we use forced end because start() could have been called more times than end()
        for (int i = 0; i < child_timers.size(); i++) {
            child_timers.at(i)->forced_end(time);
        }
    }
}

void Timer::insert_child(Timer* child) {
    child_timers.push_back(child);
}

void Timer::print(int indent_level) {
    //dont print if the timer has an empty tag (= root node)
    if (this->tag().length() > 0) {
        xprintf(Msg, "%s :    %i    %.2f \n", this->tag().c_str(), this->call_count(), this->cumulative_time());
        for (int i = 0; i < child_timers.size(); i++) {
            child_timers.at(i)->print(indent_level + 1);
        }
    }
    else {
        for (int i = 0; i < child_timers.size(); i++) {
            child_timers.at(i)->print(0);
        }
    }

}

TimerFrame::TimerFrame(string tag) {
    this->_parent = NULL;
    this->closed = false;

    map<string, TimerFrame*>::iterator i = TimerFrame::frames()->find(tag);
    if (i == TimerFrame::frames()->end()) { //not found
        TimerFrame::frames()->insert(std::make_pair(tag, this));
    }
    else {
        this->_parent = i->second;
        i->second = this;
    }

    this->tag = tag;
    Profiler::instance()->start(tag);
}

TimerFrame::~TimerFrame() {
    close();
}

void TimerFrame::close() {

    if (closed)
        return;

    closed = true;

    map<string, TimerFrame*>::iterator i = TimerFrame::frames()->find(this->tag);
    if (i != TimerFrame::frames()->end() && i->second == this) {
        if (this->parent() != NULL) {
            i->second = this->parent();
        }
        else {
            TimerFrame::frames()->erase(tag);
        }

        Profiler::instance()->end(tag);
    }
}

void TimerFrame::endTimer(string tag) {

    map<string, TimerFrame*>::const_iterator i = _frames.find(tag);
    if (i != _frames.end()) {
        i->second->close(); //manually close (we call close because it not a good idea to call destructor explicitly)
    }
}