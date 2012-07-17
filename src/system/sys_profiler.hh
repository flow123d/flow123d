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
 * @brief Class for profiling tree nodes.
 */
class Timer {
private:
    double start_time;
    double cumul_time;
    int count;
    int start_count;
    int sub_frames;
    bool running;
    string timer_tag;
    Timer* parent_timer;
    vector<Timer*> child_timers;

    void stop(double time);

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
        sub_frames = info;
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

    ~Timer();
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
    int task_size;
    string out_dir;

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
    static void initialize(MPI_Comm communicator, string odir) {
        if (!_instance)
            _instance = new Profiler(communicator);
        _instance->out_dir = odir;
    }

    /**
     * Starts a timer with specified name. If the timer is not already created, creates a new one.
     *
     * @param tag - name of the timer to start
     */
    void start(string tag);

    /**
     * Stops a timer with specified name.
     *
     * @param tag - name of the timer to stop
     */
    void end(string tag = "");

    /**
     * Sets the size of the task. Will be written into output
     *
     * @param size - size of the task
     */
    void set_task_size(int size);

    /**
     * Sets the number of subframes (eg. iterations) in which the current Timer is divided.
     *
     * @param tag - the tag of the currently running timer. If the tag doesn't match the currently
     * running one, no subframes are set.
     * @param n_subframes - the number of subframes
     */
    void set_timer_subframes(string tag, int n_subframes);
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
 */
#ifdef DEBUG_PROFILER
#define END_TIMER(tag) TimerFrame::endTimer(tag)          // only if you want end on different place then end of function
#else
#define END_TIMER(tag)
#endif

/**
 * \def SET_TIMER_SUBFRAMES(tag, subframes)
 *
 * Sets specified amount of subframes (eg. iterations) for the given tag.
 * The specified timer tag must represent the currently active timer.
 */
#ifdef DEBUG_PROFILER
#define SET_TIMER_SUBFRAMES(tag, subframes) Profiler::instance->setTimerSubframes(tag, info)
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
 */
class TimerFrame {
private:
    string tag;
    TimerFrame* _parent;
    bool closed;
    static map<string, TimerFrame*> _frames;
public:

    /**
     * Parent of the TimerFrame object (it is a TimerFrame object with the same tag,
     * but defined in the superior block of code or function)
     */
    TimerFrame* parent() {
        return _parent;
    }

    TimerFrame(string tag);

    ~TimerFrame();

    /**
     * If not already closed, closes the TimerFrame object.
     * Asks Profiler to end a timer with specified tag and changes the frames
     * map appropriately (if the TimerFrame object has a parent, associate hits parent
     * with the tag or if not, delete the tag from the map)
     */
    void close();

    /**
     * Stops the timer manually
     * @param tag - timer name
     */
    static void endTimer(string tag);

    /**
     * Tags with associated TimerFrame objects
     */
    static map<string, TimerFrame*>* frames() {
        return &_frames;
    }
};

#endif
