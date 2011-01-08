
#ifndef PROFILER_H
#define	PROFILER_H

#include <map>
#include <iostream>
#include <time.h>
#include <vector>

using namespace std;


/*
 * Class for profiling tree nodes. Maybe private for Profiler.
 */
class Timer {
private:
    double start_time;
    double cumul_time;
    int count;
    int start_count;
    bool running;
    string timer_tag;
    Timer* parent_timer;
    vector<Timer*> child_timers;

    void stop(double time);

public:
    Timer(string tag, Timer* parent);
    void start(double time); // time given by Profiler
    void end(double time);
    void forced_end(double time);
    void insert_child(Timer* child);

    int call_count() const {
        return count;
    }

    double cumulative_time() const {
        return cumul_time;
    }

    string tag() const {
        return timer_tag;
    }

    Timer* parent() {
        return parent_timer;
    }

    //print(int indent_level, MPI_Comm *comm);
    ~Timer();
};

class Profiler {
private:
    static Profiler* _instance;
    Timer *root;
    Timer *actual_node;
    clock_t start_clock;

    map<string, Timer*> tag_map;

    /**Gets the time in milliseconds since the program was launched
     */
    double inline get_time();

    /** Default constructor creates global profiling object.
     *  it reads start time of the whole program and creates:
     *  - tag to timing obj. map
     */
    Profiler();

    Profiler(Profiler const&);            // copy constructor is private

    Profiler& operator=(Profiler const&);  // assignment operator is private

    /**
     *  Pass thorugh the profiling tree (colective over processors)
     *  Print cumulative times average, balace (max/min), count (denote diferences)
     *  Destroy all structures.
     */
    ~Profiler();

public:

    //implements singleton pattern

    static Profiler* instance() {
        if (!_instance)
            _instance = new Profiler;

        return _instance;
    }

    /**
     * starts particular timing period:
     * - if the tag is new:
     *        make new instance of timing class and connect it to actual leaf of the profiling tree
     *        register into map
     *        if it is the first timing, set it as a root
     * - if tag exists, change actual leaf, increment call count of the timing object
     */
    void start(string tag);

    /**
     * - check if tag match tag of actual timing object, if not print warning and proceed to matchin parent, while closing all open timings
     * - else close actual timing and proceed to the parent node in profiling tree
     */
    void end(string tag = ""); // without tag we do not perform check

};

/*
 * Class for automatic closing.
 */
#define START_TIMER(tag)  TimerFrame(tag)
#define END_TIMER(tag) Profiler::instance()->end(tag)           // only if you want end on diferent place then end of function

class TimerFrame {
public:
    TimerFrame(string tag) {
        Profiler::instance()->start(tag);
    }

    ~TimerFrame() {
        Profiler::instance()->end();
    }
};

#endif