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
 * @file    sys_profiler.cc
 * @ingroup system
 * @brief   Profiler
 */

#include <fstream>
#include <iomanip>
#include <sys/param.h>
#include "Python.h"

#include "sys_profiler.hh"
#include "system/system.hh"
#include "Python.h"
#include "system/python_loader.hh"
#include <boost/format.hpp>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "petscsys.h" 

#include "system/file_path.hh"
#include "system/python_loader.hh"
#include "mpi.h"
#include "time_point.hh"

// namespace alias
namespace property_tree = boost::property_tree;

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
    
long MPI_Functions::sum(long* val, MPI_Comm comm) {
        long total = 0;
        MPI_Reduce(val, &total, 1, MPI_LONG, MPI_SUM, 0, comm);
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
    
long MPI_Functions::min(long* val, MPI_Comm comm) {
        long min = 0;
        MPI_Reduce(val, &min, 1, MPI_LONG, MPI_MIN, 0, comm);
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
    
long MPI_Functions::max(long* val, MPI_Comm comm) {
        long max = 0;
        MPI_Reduce(val, &max, 1, MPI_LONG, MPI_MAX, 0, comm);
        return max;
    }


#ifdef FLOW123D_DEBUG_PROFILER
/*********************************************************************************************
 * Implementation of class Timer
 */

Timer::Timer(const CodePoint &cp, int parent)
: start_time(TimePoint()),
  cumul_time(0.0),
  call_count(0),
  start_count(0),
  code_point_(&cp),
  full_hash_(cp.hash_),
  hash_idx_(cp.hash_idx_),
  parent_timer(parent),
  total_allocated_(0),
  total_deallocated_(0),
  max_allocated_(0),
  current_allocated_(0),
  alloc_called(0),
  dealloc_called(0),
  petsc_start_memory(-1),
  petsc_end_memory (-1),
  petsc_peak_memory(-1),
  petsc_local_peak_memory(-1),
  petsc_memory_difference(0)
{
    for(unsigned int i=0; i< max_n_childs ;i++)   child_timers[i]=-1;
}

void Timer::info() {
    cout << "  Timer: " << tag() << " " << this << endl;
    cout << "          start: " << petsc_start_memory << endl;
    cout << "          stop : " << petsc_end_memory << endl;
    cout << "          diff : " << petsc_memory_difference << " (" << petsc_end_memory - petsc_start_memory << ")" << endl;
    cout << "          peak : " << petsc_peak_memory << " (" << petsc_local_peak_memory << ")" << endl;
    cout << endl;
}


double Timer::cumulative_time() const {
    return cumul_time;
}

void Timer::accept_from_child(Timer &timer) {
    total_allocated_ += timer.total_allocated_ - total_allocated_;
    total_deallocated_ += timer.total_deallocated_ - total_deallocated_;
    alloc_called += timer.alloc_called - alloc_called;
    dealloc_called += timer.dealloc_called - dealloc_called;
    
    petsc_memory_difference += timer.petsc_memory_difference - petsc_memory_difference;
    
    // set current memory usage from child, and check new peak
    current_allocated_ = timer.current_allocated_;
    if (current_allocated_ > max_allocated_)
        max_allocated_ = current_allocated_;
        
    // maximum value from child are always greater or equal to parent frame
    max_allocated_ = timer.max_allocated_;
    petsc_peak_memory = timer.petsc_peak_memory
}

void Timer::accept_from_parent(Timer &timer) {
    // set default variable values from parent
    total_allocated_ += timer.total_allocated_;
    total_deallocated_ += timer.total_deallocated_;
    alloc_called += timer.alloc_called;
    dealloc_called += timer.dealloc_called;
    
    petsc_memory_difference += timer.petsc_memory_difference;
    
    // set current memory usage from parents, and check new peak
    max_allocated_ = timer.max_allocated_;
    current_allocated_ = timer.current_allocated_;
    if (current_allocated_ > max_allocated_)
        max_allocated_ = current_allocated_;
}

void Timer::pause() {
    // get the maximum resident set size (memory used) for the program.
    PetscMemoryGetMaximumUsage(&petsc_local_peak_memory);
    if (petsc_peak_memory < petsc_local_peak_memory)
        petsc_peak_memory = petsc_local_peak_memory;
}
void Timer::resume() {
    // tell PETSc to monitor the maximum memory usage so
    //   that PetscMemoryGetMaximumUsage() will work.
    PetscMemorySetGetMaximumUsage();
}

void Timer::start() {
    // Tell PETSc to monitor the maximum memory usage so
    //   that PetscMemoryGetMaximumUsage() will work.
    PetscMemorySetGetMaximumUsage();
    PetscMemoryGetCurrentUsage (&petsc_start_memory);
    
    if (start_count == 0) {
        start_time = TimePoint();
    }
    call_count++;
    start_count++;
}



bool Timer::stop(bool forced) {
    // get current memory usage
    PetscMemoryGetCurrentUsage (&petsc_end_memory);
    petsc_memory_difference += petsc_end_memory - petsc_start_memory;
    
    // get the maximum resident set size (memory used) for the program.
    PetscMemoryGetMaximumUsage(&petsc_local_peak_memory);
    if (petsc_peak_memory < petsc_local_peak_memory)
        petsc_peak_memory = petsc_local_peak_memory;
        
    if (forced) start_count=1;

    if (start_count == 1) {
        cumul_time += (TimePoint() - start_time);
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
    timers_.push_back( Timer(main_cp, -1) );
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
    Timer parent_timer = timers_[actual_node];
    
    int child_idx = find_child(cp);
    if (child_idx < 0) {
        // tag not present - create new timer
        child_idx=timers_.size();
        timers_.push_back( Timer(cp, actual_node) );
        timers_[actual_node].add_child(child_idx , timers_.back() );
    }
    actual_node=child_idx;
    
    // pause current timer
    parent_timer.pause();
    
    // propagate metrics from parent to child
    timers_[actual_node].accept_from_parent(parent_timer);
    timers_[actual_node].start();
    
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
    int child_timer = actual_node;
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
                actual_node = timers_[actual_node].parent_timer;
                
                if (actual_node >= 0) {
                    // propagate metrics from parent to child
                    timers_[actual_node].accept_from_child(timers_[child_timer]);
                    // resume current timer
                    timers_[actual_node].resume();
                }
                return;
            }
        }
        // node not found - do nothing
        return;
    }
    // node to close match the actual
    timers_[actual_node].stop(false);
    actual_node = timers_[actual_node].parent_timer;
    
    if (actual_node >= 0) {
        // propagate metrics from parent to child
        timers_[actual_node].accept_from_child(timers_[child_timer]);
        // resume current timer
        timers_[actual_node].resume();
    }
}



void Profiler::stop_timer(int timer_index) {
    if (timer_index <0) timer_index=actual_node;
    ASSERT_LESS( timer_index, timers_.size() );

    if (! timers_[timer_index].running() ) return;
    
    int child_timer = actual_node;
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
                
                if (actual_node >= 0) {
                    // propagate metrics from parent to child
                    timers_[actual_node].accept_from_child(timers_[child_timer]);
                    // resume current timer
                    timers_[actual_node].resume();
                }
                return;
            }
        // node not found - do nothing
        return;
    }
    // node to close match the actual
    timers_[actual_node].stop(false);
    actual_node=timers_[actual_node].parent_timer;
    
    if (actual_node >= 0) {
        // propagate metrics from parent to child
        timers_[actual_node].accept_from_child(timers_[child_timer]);
        // resume current timer
        timers_[actual_node].resume();
    }
}



void Profiler::add_calls(unsigned int n_calls) {
    timers_[actual_node].call_count += n_calls-1;
}


void Profiler::notify_malloc(const size_t size) {
    timers_[actual_node].total_allocated_ += size;
    timers_[actual_node].current_allocated_ += size;
    timers_[actual_node].alloc_called++;
        
    if (timers_[actual_node].current_allocated_ > timers_[actual_node].max_allocated_)
        timers_[actual_node].max_allocated_ = timers_[actual_node].current_allocated_;
    
}



void Profiler::notify_free(const size_t size) {
    timers_[actual_node].total_deallocated_ += size;
    timers_[actual_node].current_allocated_ -= size;
    timers_[actual_node].dealloc_called++;
}


double Profiler::get_resolution () {
    const int measurements = 100;
    double result = 0;

    // perform 100 measurements
    for (unsigned int i = 1; i < measurements; i++) {
        TimePoint t1 = TimePoint ();
        TimePoint t2 = TimePoint ();

        // double comparison should be avoided
        while ((t2 - t1) == 0) t2 = TimePoint ();
        // while ((t2.ticks - t1.ticks) == 0) t2 = TimePoint ();

        result += t2 - t1;
    }

    return (result / measurements) * 1000; // ticks to seconds to microseconds conversion
}



std::string common_prefix( std::string a, std::string b ) {
    if( a.size() > b.size() ) std::swap(a,b) ;
    return std::string( a.begin(), std::mismatch( a.begin(), a.end(), b.begin() ).first ) ;
}



template<typename ReduceFunctor>
void Profiler::add_timer_info(ReduceFunctor reduce, property_tree::ptree* holder, int timer_idx, double parent_time) {

    // get timer and check preconditions
    Timer &timer = timers_[timer_idx];
    ASSERT( timer_idx >=0, "Wrong timer index %d.\n", timer_idx);
    ASSERT( timer.parent_timer >=0 , "Inconsistent tree.\n");

    // fix path
    string filepath = timer.code_point_->file_;

    // if constant FLOW123D_SOURCE_DIR is defined, we try to erase it from beginning of each CodePoint's filepath
    #ifdef FLOW123D_SOURCE_DIR
        string common_path = common_prefix (string(FLOW123D_SOURCE_DIR), filepath);
        filepath.erase (0, common_path.size());
    #endif


    // generate node representing this timer
    // add basic information
    property_tree::ptree node;
    double cumul_time_sum;
    node.put ("tag",        (timer.tag()) );
    node.put ("file-path",  (filepath) );
    node.put ("file-line",  (timer.code_point_->line_) );
    node.put ("function",   (timer.code_point_->func_) );
    cumul_time_sum = reduce (timer, node);


    // statistical info
    if (timer_idx == 0) parent_time = cumul_time_sum;
    double percent = parent_time > 1.0e-10 ? cumul_time_sum / parent_time * 100.0 : 0.0;
    node.put<double> ("percent", 	percent);

    // write times children timers
    property_tree::ptree children;
    bool has_children = false;
    for (unsigned int i = 0; i < Timer::max_n_childs; i++) {
		if (timer.child_timers[i] > 0) {
			add_timer_info (reduce, &children, timer.child_timers[i], cumul_time_sum);
			has_children = true;
		}
    }

    // add children tag and other info if present
    if (has_children)
    	node.add_child ("children", children);

    // push itself to root ptree 'array'
	holder->push_back (std::make_pair ("", node));
}



void Profiler::update_running_timers() {
    for(int node=actual_node; node !=0; node = timers_[node].parent_timer)
        timers_[node].cumul_time += (TimePoint() - timers_[node].start_time);
    timers_[0].cumul_time += (TimePoint() - timers_[0].start_time);
}

template <class T>
void save_mpi_metric (property_tree::ptree &node, MPI_Comm comm, T * ptr, string name) {
    node.put (name+"-min", MPI_Functions::min(ptr, comm));
    node.put (name+"-max", MPI_Functions::max(ptr, comm));
    node.put (name+"-sum", MPI_Functions::sum(ptr, comm));
}

template <class T>
void save_nonmpi_metric (property_tree::ptree &node,  T * ptr, string name) {
    node.put (name+"-min", *ptr);
    node.put (name+"-max", *ptr);
    node.put (name+"-sum", *ptr);
}


#ifdef FLOW123D_HAVE_MPI
void Profiler::output(MPI_Comm comm, ostream &os) {
    int ierr, mpi_rank, mpi_size;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    ASSERT(ierr == 0, "Error in MPI test of rank.");
    MPI_Comm_size(comm, &mpi_size);

    // output header
    property_tree::ptree root, children;
    output_header (root, mpi_size);

    // recursively add all timers info
    // define lambda function which reduces timer from multiple processors
    // MPI implementation uses MPI call to reduce values
    auto reduce = [=] (Timer &timer, property_tree::ptree &node) -> double {
        int call_count = timer.call_count;
        double cumul_time = timer.cumulative_time () / 1000;
        
        long memory_allocated = (long)timer.total_allocated_;
        long memory_deallocated = (long)timer.total_deallocated_;
        long memory_peak = (long)timer.max_allocated_;
        
        int alloc_called = timer.alloc_called;
        int dealloc_called = timer.dealloc_called;
        
        long petsc_memory_difference = (long)timer.petsc_memory_difference;
        long petsc_peak_memory = (long)timer.petsc_peak_memory;
        
        save_mpi_metric<double>(node, comm, &cumul_time, "cumul-time");
        save_mpi_metric<int>(node, comm, &call_count, "call_count");
        
        save_mpi_metric<long>(node, comm, &memory_allocated, "memory-alloc");
        save_mpi_metric<long>(node, comm, &memory_deallocated, "memory-dealloc");
        save_mpi_metric<long>(node, comm, &memory_peak, "memory-peak");
        // 
        save_mpi_metric<int>(node, comm, &alloc_called, "memory-alloc-called");
        save_mpi_metric<int>(node, comm, &dealloc_called, "memory-dealloc-called");
        
        save_mpi_metric<long>(node, comm, &petsc_memory_difference, "memory-petsc-diff");
        save_mpi_metric<long>(node, comm, &petsc_peak_memory, "memory-petsc-peak");
        
        return MPI_Functions::sum(&cumul_time, comm);
    };

    add_timer_info (reduce, &children, 0, 0.0);
    root.add_child ("children", children);


    // create profiler output only once (on the first processor)
    // only active communicator should be the one with mpi_rank 0
    if (mpi_rank == 0) {
        /**
         * Flag to property_tree::write_json method
         * resulting in json human readable format (indents, newlines)
         */
        const int FLOW123D_JSON_HUMAN_READABLE = 1;
        // write result to stream
        property_tree::write_json (os, root, FLOW123D_JSON_HUMAN_READABLE);
    }
}


void Profiler::output(MPI_Comm comm) {
    //wait until profiling on all processors is finished
    MPI_Barrier(comm);
    update_running_timers();
    int ierr, mpi_rank, mpi_size;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    ASSERT(ierr == 0, "Error in MPI test of rank.");
    MPI_Comm_size(comm, &mpi_size);

    char filename[PATH_MAX];
    strftime(filename, sizeof (filename) - 1, "profiler_info_%y.%m.%d_%H-%M-%S.log.json", localtime(&start_time));
    json_filepath = FilePath(string(filename), FilePath::output_file);

    xprintf(MsgLog, "output into: %s\n", json_filepath.c_str());
    ofstream os(json_filepath.c_str());
    output(comm, os);
    os.close();
}

#endif /* FLOW123D_HAVE_MPI */

void Profiler::output(ostream &os) {
    // last update
    update_running_timers();

    // output header
    property_tree::ptree root, children;
    /**
     * Constant representing number of MPI processes
     * where there is no MPI to work with (so 1 process)
     */
    const int FLOW123D_MPI_SINGLE_PROCESS = 1;
    output_header (root, FLOW123D_MPI_SINGLE_PROCESS);


    // recursively add all timers info
    // define lambda function which reduces timer from multiple processors
    // non-MPI implementation is just dummy repetition of initial value
    auto reduce = [=] (Timer &timer, property_tree::ptree &node) -> double {
        int call_count = timer.call_count;
        double cumul_time = timer.cumulative_time () / 1000;
        
        long memory_allocated = (long)timer.total_allocated_;
        long memory_deallocated = (long)timer.total_deallocated_;
        long memory_peak = (long)timer.max_allocated_;
        
        int alloc_called = timer.alloc_called;
        int dealloc_called = timer.dealloc_called;
        
        long petsc_memory_difference = (long)timer.petsc_memory_difference;
        long petsc_peak_memory = (long)timer.petsc_peak_memory;
        
        save_nonmpi_metric<double>(node, &cumul_time, "cumul-time");
        save_nonmpi_metric<int>(node, &call_count, "call_count");
        
        save_nonmpi_metric<long>(node, &memory_allocated, "memory-alloc");
        save_nonmpi_metric<long>(node, &memory_deallocated, "memory-dealloc");
        save_nonmpi_metric<long>(node, &memory_peak, "memory-peak");
        
        save_nonmpi_metric<int>(node, &alloc_called, "memory-alloc-called");
        save_nonmpi_metric<int>(node, &dealloc_called, "memory-dealloc-called");
        
        save_nonmpi_metric<long>(node, &petsc_memory_difference, "memory-petsc-diff");
        save_nonmpi_metric<long>(node, &petsc_peak_memory, "memory-petsc-peak");
        
        return cumul_time;
    };

    add_timer_info (reduce, &children, 0, 0.0);
    root.add_child ("children", children);


    /**
     * Flag to property_tree::write_json method
     * resulting in json human readable format (indents, newlines)
     */
    const int FLOW123D_JSON_HUMAN_READABLE = 1;
    // write result to stream
    property_tree::write_json (os, root, FLOW123D_JSON_HUMAN_READABLE);
}


void Profiler::output() {
    char filename[PATH_MAX];
    strftime(filename, sizeof (filename) - 1, "profiler_info_%y.%m.%d_%H-%M-%S.log.json", localtime(&start_time));
    json_filepath = FilePath(string(filename), FilePath::output_file);

    xprintf(MsgLog, "output into: %s\n", this->json_filepath.c_str());
    ofstream os(json_filepath.c_str());
    output(os);
    os.close();
}

void Profiler::output_header (property_tree::ptree &root, int mpi_size) {
    update_running_timers();
    time_t end_time = time(NULL);

    const char format[] = "%x %X";
    char start_time_string[BUFSIZ] = {0};
    strftime(start_time_string, sizeof (start_time_string) - 1, format, localtime(&start_time));

    char end_time_string[BUFSIZ] = {0};
    strftime(end_time_string, sizeof (end_time_string) - 1, format, localtime(&end_time));

    // generate current run details

    root.put ("program-name",       flow_name_);
    root.put ("program-version",    flow_version_);
    root.put ("program-branch",     flow_branch_);
    root.put ("program-revision",   flow_revision_);
    root.put ("program-build",      flow_build_);
    root.put ("timer-resolution",   boost::format("%1.9f") % Profiler::get_resolution());
    // if constant FLOW123D_SOURCE_DIR is defined, we add this information to profiler (later purposes)
    #ifdef FLOW123D_SOURCE_DIR
    #endif

    // print some information about the task at the beginning
    root.put ("task-description",   task_description_);
    root.put ("task-size",          task_size_);

    //print some information about the task at the beginning
    root.put ("run-process-count",  mpi_size);
    root.put ("run-started-at",     start_time_string);
    root.put ("run-finished-at",    end_time_string);
}


void Profiler::transform_profiler_data (const string &output_file_suffix, const string &formatter) {
    PyObject * python_module;
    PyObject * convert_method;
    PyObject * arguments;
    PyObject * return_value;
    PyObject * tmp;
    int argument_index = 0;



    // debug info
    // cout << "Py_GetProgramFullPath: " << Py_GetProgramFullPath() << endl;
    // cout << "Py_GetPythonHome:      " << Py_GetPythonHome() << endl;
    // cout << "Py_GetExecPrefix:      " << Py_GetExecPrefix() << endl;
    // cout << "Py_GetProgramName:     " << Py_GetProgramName() << endl;
    // cout << "Py_GetPath:            " << Py_GetPath() << endl;
    // cout << "Py_GetVersion:         " << Py_GetVersion() << endl;
    // cout << "Py_GetCompiler:        " << Py_GetCompiler() << endl;


    // grab module and function by importing module profiler_formatter_module.py
    python_module = PythonLoader::load_module_by_name ("profiler.profiler_formatter_module");
    convert_method  = PythonLoader::get_callable (python_module, "convert" );

    //
    // def convert (json_location, output_file, formatter):
    //

    arguments = PyTuple_New (3);

    // set json path location as first argument
    tmp = PyString_FromString (json_filepath.c_str());
    PyTuple_SetItem (arguments, argument_index++, tmp);

    // set output path location as second argument
    tmp = PyString_FromString ((json_filepath + output_file_suffix).c_str());
    PyTuple_SetItem (arguments, argument_index++, tmp);

    // set Formatter class as third value
    tmp = PyString_FromString (formatter.c_str());
    PyTuple_SetItem (arguments, argument_index++, tmp);

    // execute method with arguments
    return_value = PyObject_CallObject (convert_method, arguments);
    //    cout << "calling python convert ('"<<json_filepath<<"', '"<<(json_filepath + output_file_suffix)<<"', '"<<formatter<<"')" << endl;


    if (PyBool_Check (return_value)) {
        // is boolean

        if (return_value == Py_True) {
            cout << "Python execution was successful" << endl;
        }else{
            cout << "Error when executing Python" << endl;
        }
    } else if (PyString_Check (return_value)) {
        // is string (holds error)

        char* error_msg = PyString_AsString (return_value);
        cout << "Error when executing Python: " << error_msg << endl;
    } else {
        cout << "Unknown result when executing Python: "<< endl;
    }
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

// Profiler* Profiler::instance() {
//     // if (Profiler::_instance != NULL)
//     //     Profiler::_instance = new Profiler();
//     return Profiler::_instance;
// }


#else // def FLOW123D_DEBUG_PROFILER

Profiler* Profiler::_instance = NULL;

void Profiler::initialize() {
    if (!_instance) _instance = new Profiler();
}

void Profiler::uninitialize() {
    Profiler * _instance = Profiler::instance();
    if (_instance)  {
        delete _instance;
        _instance = NULL;
    }
}


#endif // def FLOW123D_DEBUG_PROFILER
