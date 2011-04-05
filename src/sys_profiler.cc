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

#include "sys_profiler.hh"
#include "system.hh"
#include "system.hh"
#include <xio.h>
#include <boost/format.hpp>


Profiler* Profiler::_instance = NULL;
map<string, TimerFrame*> TimerFrame::_frames;

void pad_string(string *str, int length) {
    if (length > str->size())
        str->insert(str->size(), length - str->size(), ' ');
}

Profiler::Profiler(MPI_Comm comm) {
    F_ENTRY;

    task_size = 0;
    id = 0;
    communicator = comm;
    if (comm) {
        MPI_Comm_rank(PETSC_COMM_WORLD, &(id));
    }

    actual_node = root = new Timer("", NULL);
    root->start(0);

    start_clock = clock();
    start_time = time(NULL);
}

Profiler::~Profiler() {

    if (root) {
        root->forced_end(get_time()); //stop all running timers
    }

    //wait until profiling on all processors is finished
    MPI_Barrier(this->communicator);

    vector < vector<string>* >* timersInfo = new vector < vector<string>*>();
    add_timer_info(timersInfo, root, 0);

    //create profiler output only once (on the first processor)
    if (this->id == 0) {

        int maxTagLength = 0;
        int maxCallCountLength = 0;
        int maxTimeLength = 0;
        int maxMinMaxLength = 0;

        for (int i = 0; i < timersInfo->size(); i++) {
            for (int j = 0; j < timersInfo->at(i)->size(); j++) {
                string str = timersInfo->at(i)->at(j);
                int size = str.size();

                switch (j) {
                    case 0:
                        maxTagLength = max(maxTagLength, size);
                        break;
                    case 1:
                        maxCallCountLength = max(5, max(maxCallCountLength, size));
                        break;
                    case 2:
                        maxTimeLength = max(4, max(maxTimeLength, size));
                        break;
                    case 3:
                        maxMinMaxLength = max(7, max(maxMinMaxLength, size));
                        break;
                }
            }
        }

        int size;
        MPI_Comm_size(this->communicator, &size);

        time_t end_time = time(NULL);

        const char format[] = "%x %X";
        char startdest[BUFSIZ] = {0};
        strftime(startdest, sizeof (startdest) - 1, format, localtime(&start_time));

        char enddest[BUFSIZ] = {0};
        strftime(enddest, sizeof (enddest) - 1, format, localtime(&end_time));

        //create a file where the output will be written to
        FILE *out;
        const char fileformat[] = "profiler_%d-%m-%y_%H:%M:%S.out";
        char filename[PATH_MAX];
        strftime(filename, sizeof (filename) - 1, fileformat, localtime(&start_time));

        string full_fname=IONameHandler::get_instance()->get_output_file_name(filename);
        out = xfopen(full_fname.c_str(), "w+");

        //print some information about the task at the beginning
        xfprintf(out, "No. of processors: %i\n", size);
        xfprintf(out, "Task size: %i\n", task_size);
        xfprintf(out, "Start time: %s\n", startdest);
        xfprintf(out, "End time: %s\n", enddest);

        xfprintf(out, "----------------------------------------------------\n");

        string spaces = " ";
        pad_string(&spaces, maxTagLength);
        string calls = "calls";
        pad_string(&calls, maxCallCountLength);
        string time = "time";
        pad_string(&time, maxTimeLength);
        string minmax = "min/max";
        pad_string(&minmax, maxMinMaxLength);

        xfprintf(out, "%s %s %s %s subframes\n", spaces.c_str(), calls.c_str(), time.c_str(), minmax.c_str());

        for (int i = 0; i < timersInfo->size(); i++) {
            vector<string>* info = timersInfo->at(i);

            string tag = info->at(0);
            string calls = info->at(1);
            string time = info->at(2);
            string minMax = info->at(3);
            string subframes = "";
            if (info->size() > 4)
                subframes = info->at(4);

            pad_string(&tag, maxTagLength);
            pad_string(&calls, maxCallCountLength);
            pad_string(&time, maxTimeLength);
            pad_string(&minMax, maxMinMaxLength);

            xfprintf(out, "%s %s %s %s %s\n", tag.c_str(), calls.c_str(), time.c_str(), minMax.c_str(), subframes.c_str());
        }

        xfclose(out);

    }
}

void Profiler::setTimerSubframes(string tag, int n_subframes) {

    map<string, Timer*>::const_iterator i = tag_map.find(tag);

    if (i != tag_map.end() && i->second->tag() == tag) {
        i->second->subframes(n_subframes);
    }
    else
    {
        ASSERT(false, "Wrong timer tag while setting timer subframes");
    }
}

void Profiler::add_timer_info(vector<vector<string>*>* timersInfo, Timer* timer, int indent) {

    if (timer->tag().size() > 0) {
        int numproc;
        MPI_Comm_size(communicator, &numproc);

        int callCount = timer->call_count();
        int callCountMin = MPI_Functions::min(&callCount, communicator);
        int callCountMax = MPI_Functions::max(&callCount, communicator);

        double cumulTime = timer->cumulative_time() / 1000;
        double cumulTimeMin = MPI_Functions::min(&cumulTime, communicator);
        double cumulTimeMax = MPI_Functions::max(&cumulTime, communicator);
        double cumulTimeSum = MPI_Functions::sum(&cumulTime, communicator);

        string spaces = "";
        pad_string(&spaces, indent);

        vector<string>* info = new vector<string > ();
        info->push_back(spaces + timer->tag());
        info->push_back(boost::str(boost::format("%i%s") % callCount % (callCountMin != callCountMax ? "*" : "")));
        info->push_back(boost::str(boost::format("%.2f") % (cumulTimeSum / numproc)));
        info->push_back(boost::str(boost::format("%.2f") % (cumulTimeMax > 0.000001 ? cumulTimeMin / cumulTimeMax : 1)));
        if (timer->subframes() >= 0)
            info->push_back(boost::str(boost::format("%i") % timer->subframes()));

        timersInfo->push_back(info);
    }

    for (int i = 0; i < timer->child_timers_list()->size(); i++) {
        add_timer_info(timersInfo, timer->child_timers_list()->at(i), timer->tag().size() == 0 ? 0 : indent + 1);
    }
}

void Profiler::set_task_size(int size) {
    this->task_size = size;
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
    sub_frames = -1;
}

void Timer::start(double time) {
    F_ENTRY;

    count++;
    start_count++;
    if (!running) {

        start_time = time;
        running = true;
    }
}

bool Timer::end(double time) {
    //close the timer only if start() has been called as many times as end()

    start_count = MAX(start_count - 1, 0);
    if (start_count == 0) {
        this->stop(time);
        return true;
    }

    return false;
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
        for (int i = 0; i < child_timers.size(); i++) {
            child_timers.at(i)->forced_end(time);
        }
    }
}

void Timer::insert_child(Timer* child) {
    child_timers.push_back(child);
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
