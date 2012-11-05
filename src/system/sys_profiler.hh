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
 * $Id: profiler.hh 842 2011-01-08 17:58:15Z tomas.bambuch $
 * $Revision: 842 $
 * $LastChangedBy: tomas.bambuch $
 * $LastChangedDate: 2011-01-08 18:58:15 +0100 (So, 08 led 2011) $
 *
 * @file
 * 
 * 
 * TODO:
 * 1) START_TIMER(tag)  Profiler::instance()->start_timer(HASH(tag), tag, PLACE(__FILE__, __LINE__, __FUNCTION__))
 *    HASH is compile time, only pointer to tag and place string are stored
 *    NOTE: using C++11 constwxpr, we can have class for storing just data used to create new timer (see above),
 *    however actual creation of the timer has to be done runtime, since it depends on actually opened timer frame.
 * 2) Have Class to store data about point in the code, can be created at compile time.
 * 3) Profiler should take care about timer tree modifications. Timer Nodes just store data and provides access.
 * 4) START_GLOBAL_TIMER(tag) - this calls the start_timer, which creates local timer on the correct place in the hierarchy,
 *    further this timer is added to the list of global timers, this contains groups of timers with same tag, and
 *    collect/sum data from these timers in the report.
 * 5) Allow output even during calculation (not complete, but at least some thing)
 *    Report should conatin time of start as well as time of creation of the report or time from start of the program.
 *
 * 6) Every Timer node has small array( map or hash) with its own sub timers, may be better to have global array whith all
 *    Timers.
 *
 * 7) When generating report we has to deal with possibly different trees at every MPI process.
 *
 *
 * - optimize (map lookup, TimeFrame creation)
 *   - use compile time hashes instead of tags
 *   - have only one tag/hash map to Timer nodes, detect error when starting inappropriate tag (that is not child of current frame
 *     this has to be done only once - use static variable to hold pointer to appropriate (validated) Timer
 *   -    
 * - allow memory profiling 
 *   in our own new and xmalloc functions - register allocatied and deallocated memory to active Profiler frame.
 * 
 * Design:
 * basically we provide only one macro START_FRAME that should
 * - start local timer
 * - create local variable to automatically end the timer at least at return point or end of the block
 * - increase count
 * - set global pointer to current timer frame (we use it to set subframes and monitor memory allocations)
 *
 * - when called for the first time:
 *   - register to global Profiler class ( see if there is some frame of the same name with same parent)
 *   - keep pointer to parent timer frame
 */


#ifndef PROFILER_H
#define	PROFILER_H

#include <map>
#include <iostream>
#include <time.h>
#include <vector>
#include <string>

#include <mpi.h>
#include <cstring>
#include "system/const_hashes.h"

using namespace std;


class MPI_Functions {
public:

    static int sum(int* val, MPI_Comm comm) {
        int total = 0;
        MPI_Reduce(val, &total, 1, MPI_INT, MPI_SUM, 0, comm);
        return total;
    }

    static double sum(double* val, MPI_Comm comm) {
        double total = 0;
        MPI_Reduce(val, &total, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        return total;
    }

    static int min(int* val, MPI_Comm comm) {
        int min = 0;
        MPI_Reduce(val, &min, 1, MPI_INT, MPI_MIN, 0, comm);
        return min;
    }

    static double min(double* val, MPI_Comm comm) {
        double min = 0;
        MPI_Reduce(val, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
        return min;
    }

    static int max(int* val, MPI_Comm comm) {
        int max = 0;
        MPI_Reduce(val, &max, 1, MPI_INT, MPI_MAX, 0, comm);
        return max;
    }

    static double max(double* val, MPI_Comm comm) {
        double max = 0;
        MPI_Reduce(val, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        return max;
    }
};


/**
 * @brief Class that represents point in the code.
 */
/*
class CodePoint {
public:
    constexpr CodePoint(const char *tag, const char * code_point)
    : tag_(tag), code_point_(code_point), hash_(CONSTHASH(tag))
    {}

private:

    const char *tag_;
    const char *code_point_;
    unsigned int hash_;
};
*/

/**
 * @brief Class for profiling tree nodes.
 *
 * One Timer represents one particular time frame in the execution tree.
 * It collect information about total time, number of calls, allocated and deallocated memory.
 *
 * Subframes can be used to deal with calling large number of sub frames, that can not be measured itself.
 * ... this has to be solved by particular macro, that creates child timer with specified number of calls, and
 * unknown time.
 */
class Timer {

public:
    Timer(string tag, Timer* parent);
    void start(double time); // time given by Profiler
    bool end(double time);
    void forced_end(double time);
    void insert_child(Timer* child);

    /**
     * Indicates how many times the timer has been started
     */
    int call_count() const {
        return count;
    }

    int subframes() const {
        return sub_frames;
    }

    void subframes(int info) {
        sub_frames += info;
    }

    /**
     * Total time measured by this timer
     */
    double cumulative_time() const {
        return cumul_time;
    }

    /**
     * Name of the timer
     */
    string tag() const {
        return timer_tag;
    }

    /**
     * Parent of the timer
     */
    Timer* parent() {
        return parent_timer;
    }
    vector<Timer*>* child_timers_list() {
        return &child_timers;
    }


    void add_to_total_allocated(const size_t size) {
        total_allocated_ += size;
    }


    void add_to_total_deallocated(const size_t size) {
        total_deallocated_ += size;
    }

    ~Timer();

private:
    /**
     *   Start time when frame opens.
     */
    double start_time;
    /**
     * Cumulative time spent in the frame.
     */
    double cumul_time;
    /**
     * Total number of opening of the frame.
     */
    int count;
    /**
     * Number of recursive openings. ?? Would we allow this?
     */
    int start_count;
    /**
     * TODO: replace subframes with proper childs.
     */
    int sub_frames;
    /**
     * TODO: use only start_count
     */
    bool running;
    /**
     * Tag of the Timer frame. Used to identifie the frame in code and in final profiler info table.
     * Possibly we should also store compile-time computed hashes.
     */
    string timer_tag;
    /**
     * Parent in the tree.
     */
    Timer* parent_timer;
    /**
     * Should be local map or hash of child timers.
     */
    vector<Timer*> child_timers;

    /**
     * Total number of bytes allocated directly in this frame (not include subframes).
     */
    size_t total_allocated_;
    /**
     * Total number of bytes deallocated directly in this frame (not include subframes).
     */
    size_t total_deallocated_;

    void stop(double time);

};



/**
 *
 * @brief Main class for profiling by measuring time intervals.
 *
 * These time intervals form a tree structure where each interval is represented 
 * by a Timer object. The root node of the tree is automatically created and
 * started after creating the Profiler object and cannot be stopped manually.
 *
 * The class implements a singleton pattern and all the functions are accessible trough
 * Profiler::instance(), but in most cases the programmer will access the profiler
 * functions via the #START_TIMER and #END_TIMER macros. The #START_TIMER macro
 * is responsible for the fact that we don't have to call #END_TIMER macro to stop the timer and
 * the timer will be stopped at the end of the block in which #START_TIMER was used.
 * These macros internally use the TimerFrame objects and the programmer should
 * not use the TimerFrame objects directly.
 *
 * By using #SET_TIMER_SUBFRAMES macro, the programmer can specify the number of subframes (eg. iterations)
 * for the currently active timer.
 *
 */
class Profiler {
private:
    static Profiler* _instance;
    Timer *root;
    Timer *actual_node;
    clock_t start_clock;
    time_t start_time;
    MPI_Comm communicator;
    int id;

    // header informations
    int task_size_;
    string task_description_;

    string flow_name_;
    string flow_version_;
    string flow_branch_;
    string flow_revision_;
    string flow_build_;

    map<string, Timer*> tag_map;

    /**
     * Gets the time in milliseconds since the program was launched
     * @return time in milliseconds
     */
    double inline get_time();

    void add_timer_info(vector<vector<string>*>* timersInfo, Timer* timer, int indent);

    Profiler(MPI_Comm comm); // private constructor

    Profiler(Profiler const&); // copy constructor is private

    Profiler & operator=(Profiler const&); // assignment operator is private

    /**
     * Stop all timers, synchronize all processes, collect
     * profiling informations and write it to given stream.
     *
     *  Pass through the profiling tree (collective over processors)
     *  Print cumulative times average, balance (max/min), count (denote differences)
     *
     */
    void output(ostream &os);

    /**
     *  Calls output.
     *  Destroy all structures.
     */
    ~Profiler();

public:

    /**
     * Gets the Profiler object
     */
    static Profiler* instance() {
        //singleton pattern implementation
        if (!_instance)
            _instance = new Profiler(MPI_COMM_WORLD);

        return _instance;
    }

    /**
     * Destroys the Profiler object and causes that the statistics will be written to output
     */
    static void uninitialize() {
        if (_instance) {
            delete _instance;
            _instance = NULL;
        }
    }

    /**
     * Initializes the Profiler with specific MPI communicator object
     */
    static void initialize(MPI_Comm communicator) {
        if (!_instance)
            _instance = new Profiler(communicator);
    }

    /**
     * Starts a timer with specified name. If the timer is not already created, creates a new one.
     *
     * starts particular timing period:
     * - if the tag is new:
     *        make new instance of timing class and connect it to actual leaf of the profiling tree
     *        register into map
     *        if it is the first timing, set it as a root
     * - if tag exists, change actual leaf, increment call count of the timing object
     *
     *
     *
     * @param tag - name of the timer to start
     */
    Timer * start_timer(const string &tag);

    /**
     * Stops a timer with specified name.
     *
     * @param tag - name of the timer to stop
     */
    void stop_timer(const string &tag = "");

    /**
     * Sets task specific information. The string @p description with textual description of the task and the
     * number of elements of the mesh (parameter @p size). This is used for weak scaling graphs so it should
     * measure size of the task of the same type (same description).
     *
     */
    void set_task_info(string description, int size);

    /**
     * Sets informations about program version. This consists of @p program_version (includes program name), @p branch in the repository or rather full URL of the branch,
     * and SVN @p revision (or hash for GIT).
     *
     */
    void set_program_info(string program_name, string program_version, string branch, string revision, string build);


    /**
     * Sets the number of subframes (eg. iterations) in which the current Timer is divided.
     *
     * @param tag - the tag of the currently running timer. If the tag doesn't match the currently
     * running one, no subframes are set.
     * @param n_subframes - the number of subframes
     */
    void set_timer_subframes(string tag, int n_subframes);

    /**
     * Notification about allocation of given size.
     * Increase total allocated memory in current profiler frame.
     */
    void notify_malloc(const size_t size );

    /**
     * Notification about freeing memory of given size.
     * Increase total deallocated memory in current profiler frame.
     */
    void notify_free(const size_t size );
};



// These helper macros are necessary due to use of _LINE_ variable in START_TIMER macro.
#define _PASTE(a,b) a ## b
#define PASTE(a,b) _PASTE(a, b)

/**
 * \def START_TIMER(tag)
 *
 * Starts a timer with specified tag.
 *
 * In fact it creates an object named 'timer_' followed by the number of the line
 * where it has been used. This is done to avoid variable name conflicts when
 * using the macro more than once in one block of code.
 */
#ifdef DEBUG_PROFILER
#define START_TIMER(tag) TimerFrame PASTE(timer_,__LINE__) = TimerFrame(tag)
#else
#define START_TIMER(tag)
#endif

/**
 * \def END_TIMER(tag)
 *
 * Ends a timer with specified tag.
 * Use only if you want end on different place then end of function
 */
#ifdef DEBUG_PROFILER
#define END_TIMER(tag) Profiler::instance()->stop_timer(tag)
#else
#define END_TIMER(tag)
#endif

/**
 * \def END_START_TIMER(tag)
 *
 * Ends current timer and starts the new with given tag.
 */
#ifdef DEBUG_PROFILER
#define END_START_TIMER(tag) Profiler::instance()->stop_timer(); TimerFrame PASTE(timer_,__LINE__) = TimerFrame(tag)
#else
#define END_START_TIMER(tag)
#endif


/**
 * \def SET_TIMER_SUBFRAMES(tag, subframes)
 *
 * Sets specified amount of subframes (eg. iterations) for the given tag.
 * The specified timer tag must represent the currently active timer.
 */
#ifdef DEBUG_PROFILER
#define SET_TIMER_SUBFRAMES(tag, subframes) Profiler::instance()->set_timer_subframes(tag, subframes)
#else
#define SET_TIMER_SUBFRAMES(tag,subfarmes)
#endif

/**
 *
 * @brief Class for automatic timer closing. This class is used by #START_TIMER macro
 * and is responsible for the fact that we don't have to call #END_TIMER macro to stop the timer,
 * the timer will be stopped at the end of the block in which #START_TIMER was used.
 * 
 * The main idea of the approach described is that the TimerFrame variable will be destroyed
 * at the end of the block where #START_TIMER macro was used. In order to work properly
 * in situations where #END_TIMER was used to stop the timer manually before (but there is still the
 * variable which will be later destroyed), we have to store references to these variables and
 * destroy them on-demand.
 *
 * TODO:
 * Should only contain pointer to the Timer. And destructor, that close the timer.
 */
class TimerFrame {
private:
    Timer* timer_handle_;
public:
    TimerFrame(const string &tag) {
        timer_handle_ = Profiler::instance()->start_timer(tag);
    }

    ~TimerFrame() {
        Profiler::instance()->stop_timer();
    }
};

#endif
