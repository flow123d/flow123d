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
#include "mpi.h"
#include "timer_data.hh"
#include "system/jzon/Jzon.h"
/*
 * These should be replaced by using boost MPI interface
 */
int MPI_Functions::sum(int* val, MPI_Comm comm) {
        int total = 0;
        MPI_Reduce(val, &total, 1, MPI_INT, MPI_SUM, 0, comm);
        return total;
    }

double MPI_Functions::sum(double* val, MPI_Comm comm) {
        double total = 0;
        MPI_Reduce(val, &total, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        return total;
    }

int MPI_Functions::min(int* val, MPI_Comm comm) {
        int min = 0;
        MPI_Reduce(val, &min, 1, MPI_INT, MPI_MIN, 0, comm);
        return min;
    }

double MPI_Functions::min(double* val, MPI_Comm comm) {
        double min = 0;
        MPI_Reduce(val, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
        return min;
    }

int MPI_Functions::max(int* val, MPI_Comm comm) {
        int max = 0;
        MPI_Reduce(val, &max, 1, MPI_INT, MPI_MAX, 0, comm);
        return max;
    }

double MPI_Functions::max(double* val, MPI_Comm comm) {
        double max = 0;
        MPI_Reduce(val, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        return max;
    }


#ifdef FLOW123D_DEBUG_PROFILER
/*********************************************************************************************
 * Implementation of class Timer
 */

Timer::Timer(const CodePoint &cp, int parent)
: start_time(TimerData::get_time()),
  cumul_time(start_time - start_time),
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



double Timer::cumulative_time() const {
    return cumul_time.to_time();
}



void Timer::start() {
    if (start_count == 0) {
        TimerData::init();
        start_time = TimerData::get_time();
    }
    call_count++;
    start_count++;
}



void Timer::update() {
    cumul_time = (cumul_time + TimerData::get_time()) - start_time;
}



bool Timer::stop(bool forced) {

    if (forced) start_count=1;

    if (start_count == 1) {
        update ();
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
        ASSERT(i!=idx, "Too many children of the timer with tag '%s'\n", tag().c_str());
        idx=i;
    }
    child_timers[idx] = child_index;
}



string Timer::code_point_str() const {
    return boost::str( boost::format("%s:%d, %s()") % code_point_->file_ % code_point_->line_ % code_point_->func_ );
}

string Timer::common_path = "";

/***********************************************************************************************
 * Implementation of Profiler
 */






Profiler* Profiler::_instance = NULL;
CodePoint Profiler::null_code_point = CodePoint("__no_tag__", "__no_file__", "__no_func__", 0);



void Profiler::initialize()
{

    if (!_instance)
        _instance = new Profiler();

}


Profiler::Profiler()
: actual_node(0),
  task_size_(1),
  start_time( time(NULL) )

{
#ifdef FLOW123D_DEBUG_PROFILER
    static CONSTEXPR_ CodePoint main_cp = CODE_POINT("Whole Program");
    timers_.push_back( Timer(main_cp, 0) );
    timers_[0].start();
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
        idx = ( (unsigned int)(idx)==Timer::max_n_childs ? 0 : idx+1 );
    } while ( (unsigned int)(idx) != cp.hash_idx_ ); // passed through whole array
    return -1;
}



void Profiler::stop_timer(const CodePoint &cp) {
#ifdef FLOW123D_DEBUG
    // check that all childrens are closed
    Timer &timer=timers_[actual_node];
    for(unsigned int i=0; i < Timer::max_n_childs; i++)
        if (timer.child_timers[i] >0)
            ASSERT( ! timers_[timer.child_timers[i]].running() , "Child timer '%s' running while closing timer '%s'.\n", timers_[timer.child_timers[i]].tag().c_str(), timer.tag().c_str());
#endif
    if ( cp.hash_ != timers_[actual_node].full_hash_) {
        // timer to close is not actual - we search for it above actual
        for(unsigned int node=actual_node; node != 0; node=timers_[node].parent_timer) {
            if ( cp.hash_ == timers_[node].full_hash_) {
                // found above - close all nodes between
                for(; (unsigned int)(actual_node) != node; actual_node=timers_[actual_node].parent_timer) {
                    xprintf(Warn, "Timer to close '%s' do not match actual timer '%s'. Force closing actual.\n", cp.tag_, timers_[actual_node].tag().c_str());
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
            if ( (unsigned int)(timer_index) == node) {
                // found above - close all nodes between
                for(; (unsigned int)(actual_node) != node; actual_node=timers_[actual_node].parent_timer) {
                    xprintf(Warn, "Timer to close '%s' do not match actual timer '%s'. Force closing actual.\n", timers_[timer_index].tag().c_str(), timers_[actual_node].tag().c_str());
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






void Profiler::add_timer_info(MPI_Comm comm, Jzon::Node* holder, int timer_idx, double parent_time) {

    Timer &timer = timers_[timer_idx];

    ASSERT( timer_idx >=0, "Wrong timer index %d.\n", timer_idx);
    ASSERT( timer.parent_timer >=0 , "Inconsistent tree.\n");

    int numproc;
    MPI_Comm_size(comm, &numproc);

    int call_count = timer.call_count;
    int call_count_min = MPI_Functions::min(&call_count, comm);
    int call_count_max = MPI_Functions::max(&call_count, comm);
    int call_count_sum = MPI_Functions::sum(&call_count, comm);

    double cumul_time = timer.cumulative_time() / 1000; // in seconds
    double cumul_time_min = MPI_Functions::min(&cumul_time, comm);
    double cumul_time_max = MPI_Functions::max(&cumul_time, comm);
    double cumul_time_sum = MPI_Functions::sum(&cumul_time, comm);

    if (timer_idx == 0) parent_time = cumul_time_sum;

    double percent = parent_time > 1.0e-10 ? cumul_time_sum / parent_time * 100.0 : 0.0;
    string code_point = timer.code_point_str();
    code_point.erase(0, Timer::common_path.size());

    // write current timer info
    Jzon::Node node = Jzon::object();
    node.add ("tag", 		(timer.tag()) );
    node.add ("code_point", (code_point) );
	node.add ("percent", 	(boost::str(boost::format("[%.1f] ") % percent)) );
    node.add ("calls", 		(boost::str(boost::format("%i%s") % call_count % (call_count_min != call_count_max ? "*" : " "))) );
    node.add ("Tmax", 		(boost::str(boost::format("%.6f") % (cumul_time_max))) );
    node.add ("max/min", 	(boost::str(boost::format("%.6f") % (cumul_time_min > 1.0e-10 ? cumul_time_max / cumul_time_min : 1))) );
    node.add ("T/calls", 	(boost::str(boost::format("%.6f") % (cumul_time_sum / call_count_sum))) );
    node.add ("Ttotal", 	(boost::str(boost::format("%.6f") % (cumul_time_sum))) );

    // write times children timers
    Jzon::Node children = Jzon::array();
    bool has_children = false;
    for (unsigned int i = 0; i < Timer::max_n_childs; i++) {
		if (timer.child_timers[i] > 0) {
			add_timer_info(comm, &children, timer.child_timers[i], cumul_time_sum);
			has_children = true;
		}
    }

    // add chilren tag and other info if present
    if (has_children)
    	node.add("children", children);

    // push itself to root Jzon::Node 'array'
	holder->add (node);
}



void Profiler::update_running_timers() {
    for(int node=actual_node; node !=0; node = timers_[node].parent_timer)
        timers_[node].update();
    timers_[0].update();
}

std::string common_prefix( std::string a, std::string b ) {
    if( a.size() > b.size() ) std::swap(a,b) ;
    return std::string( a.begin(), std::mismatch( a.begin(), a.end(), b.begin() ).first ) ;
}


void Profiler::output(MPI_Comm comm, ostream &os) {


    //wait until profiling on all processors is finished
    MPI_Barrier(comm);
    update_running_timers();

    int ierr, mpi_rank;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    ASSERT(ierr == 0, "Error in MPI test of rank.");

    int mpi_size;
    MPI_Comm_size(comm, &mpi_size);

    time_t end_time = time(NULL);

    const char format[] = "%x %X";
    char start_time_string[BUFSIZ] = {0};
    strftime(start_time_string, sizeof (start_time_string) - 1, format, localtime(&start_time));

    char end_time_string[BUFSIZ] = {0};
    strftime(end_time_string, sizeof (end_time_string) - 1, format, localtime(&end_time));


    // find common prefix in file paths
    string common_path = timers_[0].code_point_str();
    for (unsigned int i = 1; i < timers_.size(); i++)
    	common_path = common_prefix (common_path, timers_[i].code_point_str());
    Timer::common_path = common_path;

    // generate current run details
    Jzon::Node root = Jzon::object ();
    root.add ("program-name", 		flow_name_);
    root.add ("program-version", 	flow_version_);
    root.add ("program-branch", 	flow_branch_);
    root.add ("program-revision", 	flow_revision_);
    root.add ("program-build", 		flow_build_);
    root.add ("timer_resolution", 	150);

    // print some information about the task at the beginning
    root.add ("task-description",  	task_description_);
    root.add ("task-Task size",  	task_size_);

    //print some information about the task at the beginning
    root.add ("run-process-count", 	mpi_size);
    root.add ("run-started-at", 	start_time_string);
    root.add ("run-finished-at", 	end_time_string);

    // recursively add all timers info
    Jzon::Node children = Jzon::array();
    add_timer_info(comm, &children, 0, 0.0);
    root.add("children", children);

    // write result to file
    Jzon::Writer writer;
    writer.setFormat (Jzon::StandardFormat);
    writer.writeStream (root, os);
}



void Profiler::output(MPI_Comm comm) {
	char filename[PATH_MAX];
	strftime(filename, sizeof (filename) - 1, "profiler_info_%y.%m.%d_%H-%M-%S.log", localtime(&start_time));
	string full_fname =  FilePath(string(filename), FilePath::output_file);

	xprintf(MsgLog, "output into: %s\n", full_fname.c_str());
	ofstream os(full_fname.c_str());
	output(comm, os);
	os.close();
}



void Profiler::uninitialize() {
    if (_instance) {
        ASSERT( _instance->actual_node==0 , "Forbidden to uninitialize the Profiler when actual timer is not zero (but '%s').\n",
                _instance->timers_[_instance->actual_node].tag().c_str());
        _instance->stop_timer(0);
        delete _instance;
        _instance = NULL;
    }
}

#else // def FLOW123D_DEBUG_PROFILER

Profiler* Profiler::_instance = NULL;

void Profiler::initialize() {
    if (!_instance) _instance = new Profiler();
}

void Profiler::uninitialize() {
    if (_instance)  {
        delete _instance;
        _instance = NULL;
    }
}


#endif // def FLOW123D_DEBUG_PROFILER


