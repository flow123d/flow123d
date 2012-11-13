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
 * @ingroup system
 * @brief  Profiler
 *
 */

#include <fstream>
#include <iomanip>


#include <sys/param.h>

#include "sys_profiler.hh"
#include "system/system.hh"
#include <boost/format.hpp>

#include "system/file_path.hh"

#ifdef DEBUG_PROFILER
/*********************************************************************************************
 * Implementation of class Timer
 */

Timer::Timer(const CodePoint &cp, int parent)
: start_time(0),
  cumul_time(0),
  call_count(0),
  start_count(0),
  code_point_(&cp),
  full_hash_(cp.hash_),
  hash_idx_(cp.hash_idx_),
  parent_timer(parent),
  total_allocated_(0),
  total_deallocated_(0)
{
    for(unsigned int i=0; i< max_n_childs ;i++)   child_timers[i]=-1;
}


/*
 * Standard C++ clock() function measure the time spent in calling process, not the wall time,
 * that is exactly what I want. Disadvantage is that it has small resolution about 10 ms.
 *
 * For more precise wall clock timer one can use:
 * time(&begin);
 * time(&end);
 * cout << "Time elapsed: " << difftime(end, begin) << " seconds"<< endl;
 */
Timer::TimeData Timer::get_time() {
    return clock();
}



double Timer::cumulative_time() const
{
    return 1000.0 * double(cumul_time) / CLOCKS_PER_SEC;
}



void Timer::start() {
    if (start_count == 0) start_time = get_time();
    call_count++;
    start_count++;
}



void Timer::update() {
    TimeData time = get_time();
    cumul_time += time - start_time;
}



bool Timer::stop(bool forced) {

    if (forced) start_count=1;

    if (start_count == 1) {
        TimeData time = get_time();
        cumul_time += time - start_time;
        start_count--;
        return true;
    } else {
        start_count--;
    }
    return false;
}



void Timer::add_child(int child_index, const Timer &child)
{
    unsigned int idx = child.hash_idx_;
    if (child_timers[idx] >0) {
        unsigned int i=idx;
        do {
            i=( i < max_n_childs ? i+1 : 0);
        } while (i!=idx && child_timers[i] >0);
        ASSERT(i!=idx, "Too many children of the timer with tag '%s'\n", tag() );
        idx=i;
    }
    child_timers[idx] = child_index;
}




/***********************************************************************************************
 * Implementation of Profiler
 */






Profiler* Profiler::_instance = NULL;
CodePoint Profiler::null_code_point = CodePoint("__no_tag__", "__no_file__", "__no_func__", 0);



void Profiler::initialize(MPI_Comm communicator)
{

    if (!_instance)
        _instance = new Profiler(communicator);
    else
        xprintf(Warn, "The profiler already initialized.\n");

}


Profiler::Profiler(MPI_Comm comm)
: actual_node(0),
  communicator_(comm),
  task_size_(1),
  start_time( time(NULL) )

{
#ifdef DEBUG_PROFILER
    MPI_Comm_rank(communicator_, &(mpi_rank_));

    static CONSTEXPR_ CodePoint main_cp = CODE_POINT("Whole Program");
    timers_.push_back( Timer(main_cp, 0) );
#endif
}




void Profiler::set_task_info(string description, int size) {
    task_description_ = description;
    task_size_ = size;
}



void Profiler::set_program_info(string program_name, string program_version, string branch, string revision, string build) {
    flow_name_ = program_name;
    flow_version_ = program_version;
    flow_branch_ = branch;
    flow_revision_ = revision;
    flow_build_ = build;
}



int  Profiler::start_timer(const CodePoint &cp) {

    int child_idx = find_child(cp);
    if (child_idx < 0) {
        // tag not present - create new timer
        child_idx=timers_.size();
        timers_.push_back( Timer(cp, actual_node) );
        timers_[actual_node].add_child(child_idx , timers_.back() );
    }

    timers_[child_idx].start();
    actual_node=child_idx;

    return actual_node;
}



int Profiler::find_child(const CodePoint &cp) {
    Timer &timer =timers_[actual_node];
    int idx = cp.hash_idx_;
    int child_idx;
    do {
        child_idx=timer.child_timers[idx];

        if (child_idx < 0) break; // tag is not there

        ASSERT_LESS( child_idx, timers_.size());
        if (timers_[child_idx].full_hash_ == cp.hash_) return child_idx;
        idx = ( idx==Timer::max_n_childs ? 0 : idx+1 );
    } while (idx != cp.hash_idx_); // passed through whole array
    return -1;
}



void Profiler::stop_timer(const CodePoint &cp) {
#ifdef DEBUG
    // check that all childrens are closed
    Timer &timer=timers_[actual_node];
    for(unsigned int i=0; i < Timer::max_n_childs; i++)
        if (timer.child_timers[i] >0)
            ASSERT( ! timers_[timer.child_timers[i]].running() , "Child timer '%s' running while closing timer '%s'.\n", timers_[timer.child_timers[i]].tag(), timer.tag() );
#endif

    if ( cp.hash_ != timers_[actual_node].full_hash_) {
        DBGMSG("close '%s' actual '%s'\n", cp.tag_, timers_[actual_node].tag());
        // timer to close is not actual - we search for it above actual
        for(unsigned int node=actual_node; node != 0; node=timers_[node].parent_timer) {
            DBGMSG("cmp close '%s' idx '%s'\n", cp.tag_, timers_[node].tag());
            if ( cp.hash_ == timers_[node].full_hash_) {
                // found above - close all nodes between
                for(; actual_node != node; actual_node=timers_[actual_node].parent_timer) {
                    xprintf(Warn, "Timer to close '%s' do not match actual timer '%s'. Force closing actual.\n", cp.tag_, timers_[actual_node].tag());
                    timers_[actual_node].stop(true);
                }
                // close 'node' itself
                timers_[actual_node].stop(false);
                actual_node=timers_[actual_node].parent_timer;
                return;
            }
        }
        // node not found - do nothing
        return;
    }
    // node to close match the actual
    timers_[actual_node].stop(false);
    actual_node=timers_[actual_node].parent_timer;
}



void Profiler::stop_timer(int timer_index) {
    if (timer_index <0) timer_index=actual_node;
    ASSERT_LESS( timer_index, timers_.size() );

    if (! timers_[timer_index].running() ) return;

    if ( timer_index != actual_node ) {
        // timer to close is not actual - we search for it above actual
        for(unsigned int node=actual_node; node != 0; node=timers_[node].parent_timer)
            if ( timer_index == node) {
                // found above - close all nodes between
                for(; actual_node != node; actual_node=timers_[actual_node].parent_timer) {
                    xprintf(Warn, "Timer to close '%s' do not match actual timer '%s'. Force closing actual.\n", timers_[timer_index].tag(), timers_[actual_node].tag());
                    timers_[actual_node].stop(true);
                }
                // close 'node' itself
                timers_[actual_node].stop(false);
                actual_node=timers_[actual_node].parent_timer;
                return;
            }
        // node not found - do nothing
        return;
    }

    // node to close match the actual
    timers_[actual_node].stop(false);
    actual_node=timers_[actual_node].parent_timer;
}



void Profiler::add_calls(unsigned int n_calls) {
    timers_[actual_node].call_count += n_calls-1;
}



void Profiler::notify_malloc(const size_t size) {
    timers_[actual_node].total_allocated_ += size;
}



void Profiler::notify_free(const size_t size) {
    timers_[actual_node].total_deallocated_ += size;
}




void pad_string(string *str, int length) {
    if (length > str->size())
        str->insert(str->size(), length - str->size(), ' ');
}


void Profiler::add_timer_info(vector<vector<string>*>* timersInfo, int timer_idx, int indent) {

    Timer &timer = timers_[timer_idx];

    ASSERT( timer_idx >=0, "Wrong timer index %d.\n", timer_idx);
    ASSERT( timer.parent_timer >=0 , "Inconsistent tree.\n");

    int numproc;
    MPI_Comm_size(communicator_, &numproc);

    int callCount = timer.call_count;
    int callCountMin = MPI_Functions::min(&callCount, communicator_);
    int callCountMax = MPI_Functions::max(&callCount, communicator_);

    double cumulTime = timer.cumulative_time() / 1000;
    double cumulTimeMin = MPI_Functions::min(&cumulTime, communicator_);
    double cumulTimeMax = MPI_Functions::max(&cumulTime, communicator_);
    double cumulTimeSum = MPI_Functions::sum(&cumulTime, communicator_);

    string spaces = "";
    pad_string(&spaces, indent);

    vector<string>* info = new vector<string > ();
    info->push_back(spaces + string(timer.tag()) );
    info->push_back(boost::str(boost::format("%i%s") % callCount % (callCountMin != callCountMax ? "*" : "")));
    info->push_back(boost::str(boost::format("%.2f") % (cumulTimeSum / numproc)));
    info->push_back(boost::str(boost::format("%.2f") % (cumulTimeMax > 1.0e-10 ? cumulTimeMin / cumulTimeMax : 1)));


    timersInfo->push_back(info);

    for (int i = 0; i < Timer::max_n_childs; i++)
        if (timer.child_timers[i] > 0)
            add_timer_info(timersInfo, timer.child_timers[i], indent + 1);

}






void Profiler::update_running_timers() {
    for(int node=actual_node; node !=0; node = timers_[node].parent_timer)
        timers_[node].update();
    timers_[0].update();
}



void Profiler::output(ostream &os) {

    const int column_space = 3;



    //wait until profiling on all processors is finished
    MPI_Barrier(this->communicator_);

    update_running_timers();

    vector < vector<string>* >* timersInfo = new vector < vector<string>*>();
    add_timer_info(timersInfo, 0, 0);

    //create profiler output only once (on the first processor)
    if (mpi_rank_ == 0) {

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

        int mpi_size;
        MPI_Comm_size(this->communicator_, &mpi_size);

        time_t end_time = time(NULL);

        const char format[] = "%x %X";
        char start_time_string[BUFSIZ] = {0};
        strftime(start_time_string, sizeof (start_time_string) - 1, format, localtime(&start_time));

        char end_time_string[BUFSIZ] = {0};
        strftime(end_time_string, sizeof (end_time_string) - 1, format, localtime(&end_time));

        //create a file where the output will be written to

        os << "Program name: " << flow_name_ << endl
           << "Program version: " << flow_version_ << endl
           << "Program branch: " << flow_branch_ << endl
           << "Program revision: " << flow_revision_ << endl
           << "Program build: " << flow_build_ << endl << endl;


        os << "Task description: " << task_description_ << endl
           << "Task size: " << task_size_ << endl << endl;

        //print some information about the task at the beginning
        os << "Run processes count: " << mpi_size << endl;
        os << "Run started at: " << start_time_string << endl;
        os << "Run finished at: " << end_time_string << endl;

        os << setfill ('-') << setw (40) << "" << endl;
        os.fill(' ');

        // header
        os << left << setw(maxTagLength) << "tag tree" << setw(column_space) << ""
           << setw(maxCallCountLength) << "calls" << setw(column_space) << ""
           << setw(maxTimeLength) << "time" << setw(column_space) << ""
           << setw(maxMinMaxLength) << "min/max" << setw(column_space) << ""
           << "subframes" << endl;

        os << setfill ('-') << setw (40) << "" << endl;
        os.fill(' ');

        for (int i = 0; i < timersInfo->size(); i++) {
            vector<string>* info = timersInfo->at(i);

            os << left << setw(maxTagLength) << info->at(0)         << setw(column_space) << ""        // tag
               << setw(maxCallCountLength) << info->at(1)   << setw(column_space) << ""       // calls
               << setw(maxTimeLength) << info->at(2)        << setw(column_space) << ""       // time
               << setw(maxMinMaxLength) << info->at(3)      << setw(column_space) << "";     // min/max
            if (info->size() > 4)
               os << info->at(4);                           // subframes
            os << endl;
        }
    }
}



void Profiler::output() {
            char filename[PATH_MAX];
            strftime(filename, sizeof (filename) - 1, "profiler_info_%y.%m.%d_%H:%M:%S.log", localtime(&start_time));
            string full_fname =  FilePath(string(filename), FilePath::output_file);

            DBGMSG("output into: %s\n", full_fname.c_str());
            ofstream os(full_fname.c_str());
            output(os);
            os.close();
}



void Profiler::uninitialize()
{
    if (_instance) {
        ASSERT( _instance->actual_node==0 , "Forbidden to uninitialize the Profiler when actual timer is not zero (but '%s').\n",
                _instance->timers_[_instance->actual_node].tag());
        _instance->stop_timer(0);
        delete _instance;
        _instance = NULL;
    }
}

#else // def DEBUG_PROFILER

Profiler* Profiler::_instance = NULL;

void Profiler::initialize(MPI_Comm communicator) {
    if (!_instance) _instance = new Profiler();
}

void Profiler::uninitialize() {
    if (_instance)  {
        delete _instance;
        _instance = NULL;
    }
}


#endif // def DEBUG_PROFILER


