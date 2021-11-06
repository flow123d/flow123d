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
 * @file    sys_profiler.hh
 * @brief   
 * @todo
 * - START_GLOBAL_TIMER(tag) - this calls the start_timer, which creates local timer on the correct place in the hierarchy,
 *    further this timer is added to the list of global timers, this contains groups of timers with same tag, and
 *    collect/sum data from these timers in the report.
 *
 * - Allow output even during calculation (not complete, but at least some thing)
 *    Report should contain time of start as well as time of creation of the report or time from start of the program.
 *
 * - When generating report we has to deal with possibly different trees at every MPI process.
 *
 * - test memory profiling
 *   in our own new and xmalloc functions - register allocatied and deallocated memory to active Profiler frame.
 *
 * - test in parallel
 * - extended output:
 *      cas na jedno volani (jina redukce nez pro kumulativni cas, pokud je pocet volani ruzny)
 *      procenta vuci predkovi
 *      code point (az nekde na konci radky)
 *
 *
 *  !!! Unfortunately using constexpr is worse (without optimization).
 *  This is probably due to use of static variable for
 *  CodePoint, the access could be slow, and computation of hash is done only once. Actually timing results
 *  are:
 *
 *  OPTIONS     OVERHEAD (compared to call 2x clock())
 *  -g, no c++11 : 18%
 *  -g,    c++11 : 60%
 *  -O3,no c++11 : 6%
 *  -O3,   c++11 : 6%
 */

#ifndef PROFILER_H
#define	PROFILER_H

#include "global_defs.h"

#include <mpi.h>
#include <ostream>
#include <unordered_map>

namespace boost { template <class T> struct hash; }
#include <boost/functional/hash/hash.hpp>      // for hash
#include <boost/ref.hpp>
#include <boost/tuple/detail/tuple_basic.hpp>  // for get

#include <nlohmann/json.hpp>

#include "time_point.hh"
#include "petscsys.h" 
#include "simple_allocator.hh"

//instead of #include "mpi.h"
//mpi declarations follows:
class MPI_Functions {
public:
    static int sum(int* val, MPI_Comm comm);
    static double sum(double* val, MPI_Comm comm);
    static long sum(long* val, MPI_Comm comm);
    
    static int min(int* val, MPI_Comm comm);
    static double min(double* val, MPI_Comm comm);
    static long min(long* val, MPI_Comm comm);
    
    static int max(int* val, MPI_Comm comm);
    static double max(double* val, MPI_Comm comm);
    static long max(long* val, MPI_Comm comm);
};

// Assuming all compilers support constexpr
#define CONSTEXPR_ constexpr

using namespace std;


// These helper macros are necessary due to use of _LINE_ variable in START_TIMER macro.
#define _PASTE(a,b) a ## b
#define PASTE(a,b) _PASTE(a, b)



/**
 * \def START_TIMER(tag)
 *
 * @brief Starts a timer with specified tag.
 *
 * In fact it creates an static constant expression that identifies the point in the code and
 * contains tag of the involved timer and its hash. Then it creates local variable that
 * calls @p Profiler::start_timer() in constructor and @p Profiler::stop_timer() in destructor.
 * This way the timer is automatically closed at the end of current block.
 *
 * ATTENTION: This macro expands to two statements so following code is illegal:
 * @code
 *      if (some_condition) START_TIMER(tag);
 * @endcode
 */
#ifdef FLOW123D_DEBUG_PROFILER
#define START_TIMER(tag) static CONSTEXPR_ CodePoint PASTE(cp_,__LINE__) = CODE_POINT(tag); TimerFrame PASTE(timer_,__LINE__) = TimerFrame( PASTE(cp_,__LINE__) )
#else
#define START_TIMER(tag)
#endif

/**
 * \def START_TIMER_EXT (tag, subtag)
 *
 * @brief Starts a timer with specified tag and subtag.
 *
 * In fact it creates an static constant expression that identifies the point in the code and
 * contains tag and subtag of the involved timer and its hash. Then it creates local variable that
 * calls @p Profiler::start_timer() in constructor and @p Profiler::stop_timer() in destructor.
 * This way the timer is automatically closed at the end of current block.
 *
 * ATTENTION: This macro expands to two statements so following code is illegal:
 * @code
 *      if (some_condition) START_TIMER_EXT(tag, subtag);
 * @endcode
 */
#ifdef FLOW123D_DEBUG_PROFILER
#define START_TIMER_EXT(tag, subtag) static CONSTEXPR_ CodePoint PASTE(cp_,__LINE__) = CODE_POINT_EXT(tag, subtag); TimerFrame PASTE(timer_,__LINE__) = TimerFrame( PASTE(cp_,__LINE__) )
#else
#define START_TIMER_EXT(tag, subtag)
#endif

/**
 * \def END_TIMER(tag)
 *
 * @brief Ends a timer with specified tag.
 *
 * Use only if you want to end timer before the end of block. Again this expands into two lines, see ATTENTION in previous macro.
 */
#ifdef FLOW123D_DEBUG_PROFILER
#define END_TIMER(tag) static CONSTEXPR_ CodePoint PASTE(cp_,__LINE__) = CODE_POINT(tag); Profiler::instance()->stop_timer( PASTE(cp_,__LINE__) )
#else
#define END_TIMER(tag)
#endif

/**
 * \def END_START_TIMER(tag)
 *
 * Ends current timer and starts the new one with given tag.  Again this expands into two lines, see ATTENTION in previous macro.
 */
#ifdef FLOW123D_DEBUG_PROFILER
#define END_START_TIMER(tag) Profiler::instance()->stop_timer(); START_TIMER(tag);
#else
#define END_START_TIMER(tag)
#endif


/**
 * \def ADD_CALLS(n_calls)
 *
 * @brief Increase number of calls in actual timer.
 *
 * Some time you want to measure a performance of a cycle with body that is below resolution of the Timer implementation.
 * If you know number of cycles, you can use this macro in following way:
 *
 * @code
 *  START_TIMER("cycle");
 *  unsigned int i;
 *  for(i =0; i<1000000; i++) i*i*i;
 *  ADD_CALLS(i);
 *  END_TIMER("cycle");
 * @endcode
 *
 * In the profiler report you get the total time spent in the cycle, and time per one call which will be average
 * time spent in the body of the cycle.
 */
#ifdef FLOW123D_DEBUG_PROFILER
#define ADD_CALLS(n_calls) Profiler::instance()->add_calls(n_calls)
#else
#define ADD_CALLS(n_calls)
#endif


#ifdef FLOW123D_DEBUG_PROFILER
#define CUMUL_TIMER(tag) Profiler::instance()->find_timer(tag).cumulative_time()
#else
#define CUMUL_TIMER(tag)
#endif

//////////////////////////////////////////////////////////////////////////////////////////////
#ifdef FLOW123D_DEBUG_PROFILER

/**
 * Variable which represents value when no subtag was specified in CodePoint class
 */
#define PROFILER_EMPTY_SUBTAG ""

/**
 * Variable used for default value in hash process
 */
#define PROFILER_HASH_DEFAULT 0

/**
 * @brief Function for compile-time hash computation.  (Needs C++x11 standard.)
 * Input, @p str, is constant null terminated string, result is unsigned int (usually 4 bytes).
 * Function has to be recursive, since standard requires that the body consists only from the return statement.
 *
 * SALT is hash for the empty string. Currently zero for simpler testing.
 */
inline CONSTEXPR_ unsigned int str_hash(const char * str, unsigned int default_value) {
    #define SALT 0 //0xef50e38f
    return (*str == 0 ? SALT : default_value + str_hash(str+1, PROFILER_HASH_DEFAULT) * 101 + (unsigned int)(*str) );
}

/**
 * Macro to generate constexpr CodePoint object.
 */
#define CODE_POINT(tag) CodePoint(tag, __FILE__, __func__, __LINE__)

/**
 * Macro to generate constexpr CodePoint object.
 */
#define CODE_POINT_EXT(tag, subtag) CodePoint(tag, subtag, __FILE__, __func__, __LINE__)




/**
 * @brief Class that represents point in the code.
 *
 * This class allow construction at compile time. And includes the information about the code point as well
 * as the 'tag' of the timer and cimpile-time computed hashes of this 'tag'. The @p hash_ is long one with
 * very small probability of collisions - this we use for comparison of tags. The @p hash_idx_ is the long hash modulo
 * length of the array of Timer's children, this is used for fast loop up into this array that servers as a simple hash table.
 */
class CodePoint {
public:
    CONSTEXPR_ CodePoint(const char *tag, const char * file, const char * func, const unsigned int line)
    : tag_(tag), subtag_(PROFILER_EMPTY_SUBTAG), file_(file), func_(func), line_(line),
      hash_(str_hash(tag, PROFILER_HASH_DEFAULT)),
      hash_idx_( str_hash(tag, PROFILER_HASH_DEFAULT)%max_n_timer_childs )
    {};
    CONSTEXPR_ CodePoint(const char *tag, const char *subtag, const char * file, const char * func, const unsigned int line)
    : tag_(tag), subtag_(subtag), file_(file), func_(func), line_(line),
      hash_(str_hash(subtag, str_hash(tag, PROFILER_HASH_DEFAULT))),
      hash_idx_( str_hash(subtag, str_hash(tag, PROFILER_HASH_DEFAULT))%max_n_timer_childs )
    {};

    /// Size of child arrays in timer nodes.
    static const unsigned int max_n_timer_childs=13;

    /// Tag of the code point.
    const char * const tag_;

    /// Subtag of the code point.
    const char * const subtag_;

    /// file name of the code point
    const char * const file_;

    /// file name of the code point
    const char * const func_;

    /// file name of the code point
    const unsigned int line_;

    /// Full 32-bit hash of the tag ( practically no chance of collision)
    unsigned int hash_;

    /// Hash modulo size of array of timer childs ( we have to check full hash to prevent collision)
    unsigned int hash_idx_;
};



/**
 * @brief Class for profiling tree nodes.
 *
 * One Timer represents one particular time frame in the execution tree.
 * It collects information about total time, number of calls, allocated and deallocated memory.
 *
 * It should be accessed only through Profiler, which is its friend class.
 *
 * TODO: for better performance: move copy hash_ and hash_idx_ into Timer since CodePoint are in static
 * variables, that may be slower to acces.
 *
 */
class Timer {


public:
    /// Size of array @p child_timers, the hash table containing descendants in the call tree.
    static const unsigned int max_n_childs=CodePoint::max_n_timer_childs;

    /**
     * Creates the timer node object. Should not be called directly, but through the START_TIMER macro.
     */
    Timer(const CodePoint &cp, int parent);


    /**
     * Start the timer. If it is already started, just increase number of starts (recursions) and calls.
     */
    void start();

    /**
     * If number of starts (recursions) drop back to zero, we stop the timer and add the period to the cumulative time.
     * This method do not take care of its childs (it has no access to the other timers).
     * When the parameter 2p forced is 'true', we stop the timer immediately regardless the number of recursions.
     * Returns true if the timer is not closed (recursions didn't drop to zero yet).
     */
    bool stop(bool forced = false);


    /// Getter for the 'tag'.
    inline string tag() const {
        string buf(code_point_->tag_);
        buf.append(code_point_->subtag_);
        return buf;
    }

    /// Returns true if the timer is open, number of starts (recursions) is nonzero.
    inline bool running() const
        { return start_count >0; }

    /// Returns string with description of the code point where the timer was first started.
    std::string code_point_str() const;

    /**
     * Returns cumulative time of the timer in seconds.
     */
    double cumulative_time() const;

    /*
     * Adds given index @p child_index of the timer @p child to the correct place in the hash table.
     */
    void add_child(int child_index, const Timer &child);


protected:
    
    /**
     * Pauses current timer, save measured petsc memory information util resume.
     * We get Petsc maximum memory usage via PetscMemoryGetMaximumUsage call
     * and save this value into temp value. (we override local maximum if temp
     * value is greater)
     */
    void pause();
    /**
     * Resume current timer. e tell Petsc to monitor the maximum memory 
     * usage once again. We call PetscMemorySetGetMaximumUsage so later in 
     * resume() method will PetscMemoryGetMaximumUsage method work.
     */
    void resume();
    
    /**
     *   Start time when frame opens.
     */
    TimePoint start_time;
    /**
     * Cumulative time spent in the frame.
     */
    double cumul_time;
    /**
     * Total number of opening of the frame.
     */
    unsigned int call_count;
    /**
     * Number of recursive openings.
     */
    unsigned int start_count;


    /**
     * Code point of the first START_TIMER for the particular tag. The 'tag' identifies timer
     * and is used in reported profiler table.
     */
    const CodePoint *code_point_;
    /// Full tag hash. Copy from code_point_.
    unsigned int full_hash_;
    /// Hash modulo size of array of timer childs. Copy from code_point_.
    unsigned int hash_idx_;

    /**
     * Index of the parent timer node  in the tree. Negative value means 'not set'.
     */
    int parent_timer;
    /**
     * Indices of the child timers in the Profiler::timers_ vector. Negative values means 'not set'.
     */
    int child_timers[max_n_childs];

    /**
     * Total number of bytes allocated in this frame. After 
     * Profiler::propagate_timers call will also contain values from children.
     */
    size_t total_allocated_;
    /**
     * Total number of bytes deallocated in this frame. After 
     * Profiler::propagate_timers call, will also contain values from children.
     */
    size_t total_deallocated_;
    /**
     * Maximum number of bytes allocated at one time in this frame. After 
     * Profiler::propagate_timers call, maximum value will be taken from this
     * Timer and also from all children Timers.
     */
    size_t max_allocated_;
    /**
     * Current number of bytes allocated in this frame at the same time. 
     * This value is used to monitor maximum bytes allocated. When notify_free
     * and notify_malloc is called this values is changed and new maximum 
     * is tested.
     */
    size_t current_allocated_;
    
    /**
     * Number of times new/new[] operator was used in this scope
     */
    int alloc_called;
    /**
     * Number of times delete/delete[] operator was used in this scope
     */
    int dealloc_called;
    
    #ifdef FLOW123D_HAVE_PETSC
    /**
     * Number of bytes used by Petsc at the start of time-frame
     */
    PetscLogDouble petsc_start_memory;
    /**
     * Number of bytes used by Petsc at the end of time-frame
     */
    PetscLogDouble petsc_end_memory;
    /**
     * Difference between start and end of a petsc memory usage 
     */
    PetscLogDouble petsc_memory_difference;
    /**
     * Maximum amount of memory used that was PetscMalloc()ed at any time 
     * during this run.
     *
     * The memory usage reported here includes all Fortran arrays (that may be
     * used in application-defined sections of code).
     */
    PetscLogDouble petsc_peak_memory;
    /**
     * Local maximum amount of memory used that was PetscMalloc()ed 
     * used during time-frame pause/resume. Auxilary variable for storing 
     * local memory used when pause is called.
     */
    PetscLogDouble petsc_local_peak_memory;
    #endif // FLOW123D_HAVE_PETSC
    
    friend class Profiler;
    friend std::ostream & operator <<(std::ostream&, const Timer&);
    
    /**
     * if under unit testing, specify friend so protected members can be tested
     */
    #ifdef __UNIT_TEST__
        friend ProfilerTest;
    #endif /* __UNIT_TEST__ */

};

/*
struct SimpleTranslator {
    typedef std::string internal_type;
    typedef int         external_type;

    // Converts a string to int
    boost::optional<external_type> get_value(const internal_type& str) {
        return boost::optional<external_type>(std::stoi(str));
    }

    // Converts a bool to string
    boost::optional<internal_type> put_value(const external_type& i){
        return boost::optional<internal_type>(std::to_string(i));
    }
};

namespace boost {
namespace property_tree {

template<typename Ch, typename Traits, typename Alloc>
struct translator_between<std::basic_string< Ch, Traits, Alloc >, int> {
    typedef SimpleTranslator type;
};


} // namespace property_tree
} // namespace boost
*/
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
 *
 * Currently the Profiler system is not thread safe. No idea how to do this.
 *
 */
class Profiler {
public:
    /**
     * Returns unique Profiler object.
     * if clear flag is set, will delete profiiler isntance
     */
    static Profiler* instance(bool clear = false);
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
     * Starts a timer with code point, tag and hashes specified by CodePoint object @p cp.
     * If the timer is not already created, it creates a new one. It returns index of
     * the actual timer.
     */
    int start_timer(const CodePoint &cp);
    /**
     * Stops actual timer. It check if the hash of given code point match hash of the
     * tag of actual timer node. If not we print out warning and try to find the correct tag
     * towards the tree root closing all nodes we pass through.
     *
     * If FLOW123D_DEBUG is set, we check that all children are closed.
     */
    void stop_timer(const CodePoint &cp);

    /**
     * Stop timer with index given by @p timer_index. If this is not equal to @p actual_node, we
     * traverse the tree towards root while force closing nodes by the way.
     *
     * Negative @p timer_index means close @p actual_node
     */
    void stop_timer(int timer_index = -1);

    /**
     * Adds @p n_calls - 1 to the total number of calls of the current timer. Minus one, since one call is counted when
     * timer was started. You should use macro ADD_CALLS above.
     */
    void add_calls(unsigned int n_calls);
    /**
     * Notification about allocation of given size.
     * Increase total allocated memory in current profiler frame.
     */
    void notify_malloc(const size_t size, const long p);
    /**
     * Notification about freeing memory of given size.
     * Increase total deallocated memory in current profiler frame.
     */
    void notify_free(const long p);

    /**
     * Return average profiler timer resolution in seconds
     * based on 100 measurements
     */
    static double get_resolution ();

    /**
     * Find a first timer matching the tag.
     * O(n) complexity.
     */
    Timer find_timer(string tag);




#ifdef FLOW123D_HAVE_MPI
    /**
     * @brief Output current timing information into the given stream.
     *
     * COLECTIVE - all processes in the communicator have to call this
     * method. All timers are finished,  all processes are synchronized, collect
     * profiling informations are collected and written to the given stream.
     *
     *  Pass through the profiling tree (collective over processors)
     *  Print cumulative times average, balance (max/min), count (denote differences)
     *
     */
    void output(MPI_Comm comm, std::ostream &os);

    /**
     * Same as previous, but output to the file with default name: "profiler_info_YYMMDD_HH::MM:SS.log".
     * Empty body if macro FLOW123D_DEBUG_PROFILER is not defined.
     */
    void output(MPI_Comm comm, string profiler_path = "");

#endif /* FLOW123D_HAVE_MPI */
    /**
     * @brief Output current timing information into the given stream.
     *
     * It temporally stops all timers, synchronize all processes, collect
     * profiling informations and write it to the given stream.
     *
     *  Pass through the profiling tree (collective over processors)
     *  Print cumulative times average, balance (max/min), count (denote differences)
     *
     */
    void output(std::ostream &os);

    /**
     * Same as previous, but output to the file with default name: "profiler_info_YYMMDD_HH::MM:SS.log".
     * Empty body if macro FLOW123D_DEBUG_PROFILER is not defined.
     */
    void output(string profiler_path = "");

    /**
     * Method will transform last profiler json file to desired format
     */
    void transform_profiler_data (const string &output_file_suffix, const string &formatter);
    /**
     * Stop all timers and destroys the Profiler object.
     * If you want some output call @p output method just before.
     */
    static void uninitialize();

    /**
     * Class-specific allocation function new. Called by the usual
     * single-object new-expressions if allocating an object of type Profiler.
     */
    static void* operator new (size_t sz);
    /**
     * Class-specific allocation function delete. Deallocates storage 
     * previously allocated by a matching operator new. These deallocation
     * functions are called by delete-expressions.
     */
    static void operator delete (void* p);
    /// Sized deallocator, doesthe same as operator delete (void* p)
    static void operator delete (void* p, std::size_t);
    
    /**
     * Public setter to turn on/off memory monitoring
     * @param global_monitor whether to turn global monitoring on or off
     * @param petsc_monitor petsc monitoring
     */
    void static set_memory_monitoring(const bool global_monitor, const bool petsc_monitor);
    
    /**
     * Public getter to memory monitoring
     * @return memory monitoring status
     */
    bool static get_global_memory_monitoring();
    
    /**
     * Public getter to petsc memory monitoring
     * @return memory monitoring status
     */
    bool static get_petsc_memory_monitoring();

    /**
     * Run calibration frame "UNIT PAYLOAD".
     * That should be about 100x timer resolution.
     */
    void calibrate();

    /**
     * Time of a unit payload, result of a single measurement. Can be used for raw calibration.
     */
    double calibration_time() {
        if (calibration_time_ < 0) calibrate();
        return calibration_time_;
    }
    /**
     * if under unit testing, specify friend so protected members can be tested
     */
    #ifdef __UNIT_TEST__
        friend ProfilerTest;
    #endif /* __UNIT_TEST__ */


protected:
    
    /**
     * Whether to monitor operator 'new/delete'
     */
    static bool global_monitor_memory;
    
    /**
     * Whether to monitor petsc memory usage
     */
    static bool petsc_monitor_memory;
    
    /**
     * When creating Profiler also reserve some bytes in malloc_map so overhead 
     * of creating single items is lowered. This value is passed as parameter in
     * map.reserve() method so it indicates how many objects (pointers) are 
     * allocated at first.
     */
    static const long malloc_map_reserve;
    

    /**
     * Method will propagate values from children timers to its parents
     */
    void propagate_timers ();
    
    /**
     * Method for exchanging metrics from child timer to its parent timer
     */
    void accept_from_child (Timer &parent, Timer &child);
    
    /**
     * Try to find timer with tag (in fact only its 32-bit hash) from given code point @p cp.
     * Returns -1 if it is not found otherwise it returns its index.
     */
    int find_child(const CodePoint &cp);


    /**
     * Method will prepare construct specific details about the run (time start and time end)
     * and write them along with basic informations about the run (name, description, ...)
     * into nlohmann:json object
     */
    void output_header (nlohmann::json &root, int mpi_size);

    /**
     * Open a new file for profiler output with default name based on the
     * actual time and date. Returns a pointer to the stream of the output file.
     */
    std::shared_ptr<std::ostream> get_default_output_stream();

    /// Vector of all timers. Whole tree is stored in this array.
    vector<Timer, internal::SimpleAllocator<Timer>> timers_;

    /// Index of the actual timer node. Negative value means 'unset'.
    unsigned int actual_node;

    /// MPI communicator used for final reduce of the timer node tree.
    //MPI_Comm communicator_;
    /// MPI_rank
    //int mpi_rank_;

    /**
     * flag indicating that collection of timer details will be
     * using MPI
    bool mpi_used;
     */
    // header informations

    /// Some measure of the size of the task in the set of the tasks that differs
    /// only by size - used for scaling tests.
    int task_size_;
    /// Task description and identifier in possible database of all Profiler results.
    string task_description_;
    /// Time and date of the start of the task solution. In fact start of the Profiler.
    time_t start_time;

    /// Name of the program.
    string flow_name_;
    /// Version of the program.
    string flow_version_;
    /// Http address of the branch in a repository.
    string flow_branch_;
    /// Revision or GIT hash.
    string flow_revision_;
    /// Build date and time.
    string flow_build_;
    /// Variable which stores last json log filepath
    string json_filepath;

    Timer none_timer_;

    /// Time of a unit payload, result of single measurement. Can be used for raw calibration.
    double calibration_time_;

protected:
    /**
     * Use DFS to pass through the tree and collect information about all timers reduced from the processes in the communicator.
     * For every timer the information strings are stored in the struct TimerInfo in order to pad fields correctly
     * to have alligned columns on the output. The alligning is performed in the output() method.
     */
    template<typename ReduceFunctor>
    void add_timer_info(ReduceFunctor reduce, nlohmann::json* node, int timer_idx, double parent_time);

    //Profiler(MPI_Comm comm); // private constructor
    Profiler(); // private constructor
    Profiler(Profiler const&); // copy constructor is private
    Profiler & operator=(Profiler const&); // assignment operator is private

};






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
    int const timer_index_;
public:
    inline TimerFrame(const CodePoint &cp)
    : timer_index_( Profiler::instance()->start_timer(cp) )
    {}

    ~TimerFrame() {
        Profiler::instance()->stop_timer(timer_index_);
    }
};


/**
 * Simple class providing static map variable storing address and alloc size
 */
// gcc version 4.9 and lower has following bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=59751
// fix in version 4.9: https://gcc.gnu.org/gcc-4.9/changes.html#cxx
// typedef unordered_map<long, int, hash<long>, equal_to<long>, internal::SimpleAllocator<pair<const long, int>>> unordered_map_with_alloc;
typedef std::unordered_map<long, int, boost::hash<long>, equal_to<long>, internal::SimpleAllocator<std::pair<const long, int>>> unordered_map_with_alloc;
class MemoryAlloc {
public:
    /**
     * Create static map containing <allocation address, allocation size> pairs
     *   map is used for storing allocations and deallocations of all object not 
     *   related to profiler after profiler initialization phase
     */
	static unordered_map_with_alloc & malloc_map();
};




#else // FLOW123D_DEBUG_PROFILER


// dummy declaration of Profiler class
class Profiler {
public:
    static Profiler* instance(bool clear = false);

    void set_task_info(string, int)
    {}
    void set_program_info(string, string, string, string, string)
    {}
    void notify_malloc(const size_t )
    {}
    void notify_free(const size_t )
    {}
    void output(MPI_Comm, ostream &)
    {}
    void output(MPI_Comm, string)
    {}
    void output(string)
    {}
    void output(MPI_Comm)
    {}
    void output()
    {}
    void transform_profiler_data(const string &, const string &)
    {}
    double get_resolution () const
    { return 0.0; }
    const char *actual_tag() const
    { return NULL; }
    inline unsigned int actual_count() const
    { return 0; }
    inline double actual_cumulative_time() const
    { return 0.0; }
    static void uninitialize();
    void calibrate();
    double calibration_time() {
        return -2;
    }
private:
    Profiler() {}
};




#endif


#endif

