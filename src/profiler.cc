
#include <sys/param.h>

#include "profiler.h"

Profiler* Profiler::_instance = NULL;

Profiler::Profiler() {
    start_clock = clock();
}

void Profiler::start(string tag) {
    map<string, Timer*>::const_iterator i = Profiler::tag_map.find(tag);

    if (i == Profiler::tag_map.end()) //the key has not been found
    {
        Timer *tmr = new Timer(tag, actual_node);
        tmr->start(get_time());
        Profiler::tag_map.insert(std::make_pair(tag, tmr));

        if (actual_node)
            actual_node->insert_child(tmr);

        if (!root)
            root = tmr;

        actual_node = tmr;
    } else {
        i->second->start(get_time());

        actual_node = i->second;
    }
}

void Profiler::end(string tag) {
    if (actual_node) {
        if (tag == "" || actual_node->tag() == tag) {
            //close actual timing and all of its children
            actual_node->end(get_time());
        } else {
            map<string, Timer*>::const_iterator i = Profiler::tag_map.find(tag);
            if (i != Profiler::tag_map.end()) {
                i->second->end(get_time());
            } else {
                //ERROR - the key has not been found
            }
        }
    } else {
        //we want to end a timing object, but there is not an active one...
    }
}

double Profiler::get_time() {
    return (double) (clock() - start_clock) / CLOCKS_PER_SEC * 1000;
}

Timer::Timer(string tag, Timer* parent) {
    this->timer_tag = tag;
    running = false;
    this->parent_timer = parent;
}

void Timer::start(double time) {
    count++;
    start_count++;
    if (!running) {
        start_time = time;
        running = true;
    }
}

void Timer::end(double time) {
    //close the timer only if start() has been called as many times as end()
    start_count = MAX(start_count - 1, 0);
    if (start_count == 0) {
        this->stop(time);
    }
}

void Timer::forced_end(double time) {
    start_count = 0;
    this->stop(time);
}

void Timer::stop(double time) {
    if (running && start_count == 0) {
        running = false;
        cumul_time += time - start_time;
        //force to close all of its children
        //we use forced end because start() could have been called more times than end()
        for (int i = 0; i < child_timers.size(); i++){
            child_timers.at(i)->forced_end(time);
        }
    }
}

void Timer::insert_child(Timer* child) {
    child_timers.push_back(child);
}
