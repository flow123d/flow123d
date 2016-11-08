#include <stdio.h> 
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <chrono>
#include <iostream>
#include <fstream>
#include <stdarg.h>
#include "libs/json.hpp"

#define KB  1024
#define MB  1048576
#define REP (MB * 3)
#define OFFSET 256
#define ARR_SIZE (MB * 100)

#define GIGA  1.0e+9
#define MEGA  1.0e+6
#define KILO  1.0e+3
#define MILI  1.0e-3
#define MICRO 1.0e-6
#define NANO  1.0e-9

#define CHAR_SIZE sizeof(char)
#define INT_SIZE  sizeof(int)

#define SHOW_DURATION true
#define SHOW_DETAILS  true

using namespace std;
using json = nlohmann::json;

class Timer {
protected:
    chrono::high_resolution_clock::time_point _start;
    chrono::high_resolution_clock::time_point _stop;
public:
    chrono::duration<double, nano> duration;
    void start() {
        this->_start = std::chrono::high_resolution_clock::now();
    }
    void stop() {
        this->_stop = std::chrono::high_resolution_clock::now();
        this->duration = chrono::duration_cast<chrono::nanoseconds>(this->_stop - this->_start);
    }
};

void printf_debug(const char * fmt, ...) {
    string newfmt = string(fmt) + "\r";
    const char * newfmt_c = newfmt.c_str();
    va_list args;
    va_start(args, fmt);
    vprintf(newfmt_c, args);
    va_end(args);
    cout << flush;
}

void cpu_test(json &results, int repetition = REP) {
    Timer timer;
    int i, sum;
    
    printf("-- running CPU REG test with %d k repetition\n", repetition / KB);
    
    //-------------------------------------------------------
    timer.start();
        for (i = 0; i < repetition; i++) {
            sum += i;
        }
    timer.stop();
    //-------------------------------------------------------
    
    if (SHOW_DURATION) {
        results["duration"] = timer.duration.count();
    }
    
    if (SHOW_DETAILS) {
        results["reps"] = repetition;
    }
    
    results["MB_per_sec"] = (repetition * INT_SIZE)/ (timer.duration.count() * NANO) / MEGA;
}


template <int N, int M>
void cpu_test_r(json &results, int (&arr)[M], int (&sizes)[N], int repetition = REP) {
    Timer timer;
    int i, j, mod, sum;
    const int i_max = sizeof(sizes) / INT_SIZE;
    
    printf("-- running CPU R   test with %d k repetition\n", repetition / KB);
    
    for (i = 0; i < i_max; i++) {
        mod = sizes[i] - 1;
        printf_debug("buffer size %9d kB...", sizes[i] / KB);
        
        //-------------------------------------------------------
        timer.start();
            for (j = 0; j < repetition; j++)  {
                sum += arr[(j * OFFSET) & mod];
            }
        timer.stop();
        //-------------------------------------------------------
        
        if (SHOW_DURATION) {
            results["duration"][to_string(sizes[i]/1024 * INT_SIZE)] = timer.duration.count() * NANO;
        }
        
        results["MB_per_sec"][to_string(sizes[i]/1024 * INT_SIZE)] = (repetition * INT_SIZE)/ (timer.duration.count() * NANO) / MEGA;
    }
    
    if (SHOW_DETAILS) {
        results["reps"] = repetition;
        results["size"] = sizeof(arr) / INT_SIZE;
        results["sum"] = sum;
    }
}


template <int N, int M>
void cpu_test_w(json &results, int (&arr)[M], int (&sizes)[N], int repetition = REP) {
    Timer timer;
    int i, j, mod;
    const int i_max = sizeof(sizes) / INT_SIZE;
    
    printf("-- running CPU W   test with %d k repetition\n", repetition / KB);
    
    for (i = 0; i < i_max; i++) {
        mod = sizes[i] - 1;
        printf_debug("buffer size %9d kB...", sizes[i] / KB);
        
        //-------------------------------------------------------
        timer.start();
            for (j = 0; j < repetition; j++) {
                arr[(j * OFFSET) & mod] = 0;
            }
        timer.stop();
        //-------------------------------------------------------
        
        if (SHOW_DURATION) {
            results["duration"][to_string(sizes[i]/1024 * INT_SIZE)] = timer.duration.count() * NANO;
        }
        
        results["MB_per_sec"][to_string(sizes[i]/1024 * INT_SIZE)] = (repetition * INT_SIZE)/ (timer.duration.count() * NANO) / MEGA;
    }
    
    if (SHOW_DETAILS) {
        results["reps"] = repetition;
        results["size"] = sizeof(arr) / INT_SIZE;
    }
}


template <int N, int M>
void cpu_test_rw(json &results, int (&arr)[M], int (&sizes)[N], int repetition = REP) {
    Timer timer;
    int i, j, mod, sum;
    const int i_max = sizeof(sizes) / INT_SIZE;
    
    printf("-- running CPU R/W test with %d k repetition\n", repetition / KB);

    for (i = 0; i < i_max; i++) {
        mod = sizes[i] - 1;
        printf_debug("buffer size %9d kB...", sizes[i] / KB);
        
        //-------------------------------------------------------
        timer.start();
            for (j = 0; j < repetition; j++) {
                arr[(j * OFFSET) & mod]++;
                arr[(j * OFFSET) & mod]--;
            }
        timer.stop();
        //-------------------------------------------------------
        
        if (SHOW_DURATION) {
            results["duration"][to_string(sizes[i]/1024 * INT_SIZE)] = timer.duration.count() * NANO;
        }
        
        results["MB_per_sec"][to_string(sizes[i]/1024 * INT_SIZE)] = (repetition * INT_SIZE)/ (timer.duration.count() * NANO) / MEGA;
    }
    
    if (SHOW_DETAILS) {
        results["reps"] = repetition;
        results["size"] = sizeof(arr) * INT_SIZE;
    }
}

template <int N>
void io_test_rw(json &results, int (&sizes)[N], int file_size, const int buffer_size = 1 * MB) {
    char *buffer = new char[buffer_size];
    char *read_buffer;
    Timer timer_r, timer_w;
    ofstream write_stream;
    ifstream read_stream;
    string fname;
    int repetition, i, j, k, i_max;
    i_max = sizeof(sizes) / INT_SIZE;
    
    printf("-- running  IO R/W bandwidth test with %1.3f MB file size\n", (float)file_size / MB);
    
    for (i = 0; i < i_max; i++) {
        printf_debug("(%2d/%2d) buffer size %9.3f kB...", i + 1, i_max, (float)sizes[i] / KB);
        repetition = file_size / (sizes[i] * CHAR_SIZE);
        fname = "tmp_file_" + to_string(sizes[i]) + ".tmp";
        read_buffer = new char[sizes[i]];
        
        // write file
        write_stream.open(fname.c_str(), ios::binary | ios::out);
        
        //-------------------------------------------------------
        timer_w.start();
            for (int j = 0; j < repetition; j++) {
                write_stream.write (buffer, sizes[i]);
            }
        timer_w.stop();
        //-------------------------------------------------------
        
        if (write_stream.fail()) {
            // buffer may be too big so ignore this 
            write_stream.close();    
            remove(fname.c_str());
            continue;
        } else {
            write_stream.close();    
        }
        
        // read file
        read_stream.open(fname.c_str(), ios::binary | ifstream::in);
        
        //-------------------------------------------------------
        timer_r.start();
            for (int j = 0; j < repetition; j++) {
                read_stream.read (read_buffer, sizes[i]);
            }
        timer_r.stop();
        //-------------------------------------------------------
        
        if (read_stream.fail()) {
            // buffer may be too big so ignore this 
            read_stream.close();    
            remove(fname.c_str());
            continue;
        } else {
            read_stream.close();    
        }
        
        remove(fname.c_str());
        
        if (SHOW_DURATION) {
            results["write"]["duration"][to_string(sizes[i] * CHAR_SIZE)] = timer_w.duration.count() * NANO;
            results["read"]["duration"][to_string(sizes[i] * CHAR_SIZE)] = timer_r.duration.count() * NANO;
        }
        
        results["write"]["MB_per_sec"][to_string(sizes[i] * CHAR_SIZE)] = ((repetition * sizes[i])/timer_w.duration.count() / NANO) / MEGA;
        results["read"]["MB_per_sec"][to_string(sizes[i] * CHAR_SIZE)] = ((repetition * sizes[i])/timer_r.duration.count() / NANO) / MEGA;
    }
    
    if (SHOW_DETAILS) {
        results["size"] = file_size;
    }
}

void io_test_many(json &results, int no_files = KILO) {
    Timer timer_r, timer_w, timer_d;
    ofstream write_stream;
    ifstream read_stream;
    string fname;
    int i, n;
    
    printf("-- running IO  R/W latency test with %d files\n", no_files);
    
    printf_debug("writing test...");
    
    //-------------------------------------------------------
    timer_w.start();
    for (i = 0; i < no_files; i++) {
        fname = "tmp_file_" + to_string(i) + ".tmp";
        write_stream.open(fname.c_str(), ios::binary | ios::out);
        write_stream << 0;
        write_stream.close();
    }
    timer_w.stop();
    //-------------------------------------------------------
    
    
    printf_debug("reading test...");
    
    //-------------------------------------------------------
    timer_r.start();
    for (i = 0; i < no_files; i++) {
        fname = "tmp_file_" + to_string(i) + ".tmp";
        read_stream.open(fname.c_str(), ios::binary | ios::in);
        read_stream >> n;
        read_stream.close();
    }
    timer_r.stop();
    //-------------------------------------------------------
    
    
    printf_debug("removing test...");
    
    //-------------------------------------------------------
    timer_d.start();
    // remove files
    for (i = 0; i < no_files; i++) {
        fname = "tmp_file_" + to_string(i) + ".tmp";
        remove(fname.c_str());
    }
    timer_d.stop();
    //-------------------------------------------------------
    
    if (SHOW_DURATION) {
        results["write"]["duration"]        = timer_w.duration.count() * NANO;
        results["read"]["duration"]         = timer_r.duration.count() * NANO;
        results["del"]["duration"]          = timer_d.duration.count() * NANO;
    }
    
    results["write"]["k_files_per_sec"]   = (float)no_files / (timer_w.duration.count() * NANO) / KILO;
    results["read"]["k_files_per_sec"]    = (float)no_files / (timer_r.duration.count() * NANO) / KILO;
    results["del"]["k_files_per_sec"]     = (float)no_files / (timer_d.duration.count() * NANO) / KILO;
    
    if (SHOW_DETAILS) {
        results["count"] = no_files;
    }
}


/**
 * Start benchmark, usage:
 * optional <output> file: will create json file output with results
 * optional <scale> float: scales repetition count for benchmark
 *                         default is 1, for example value 2 will run tests 
 *                         twice as many times, value 0.5 will experiments half
 *                         as many times
 */
int main(int argc,  char* argv[]) {
    map<int, long> results_write, results_read, results_rw, results_cpu;
    int rep_cnt = (int)(argc >= 3 ? std::stof(argv[2]) * REP : 1 * REP);
    
    // chunk size for testing
    static int sizes[] = {
        // 1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
        1 * KB, 2 * KB, 4 * KB, 8 * KB, 16 * KB, 32 * KB, 64 * KB, 128 * KB,
        256 * KB, 512 * KB, 1 * MB, 2 * MB, 4 * MB, 8 * MB, 16 * MB, 32 * MB
        // 4 * KB, 128 * KB, 8 * MB, 32 * MB
    };
    
    // chunk size for testing
    static int io_sizes[] = {
        16, 32, 64, 128, 256, 512,
        1 * KB, 2 * KB, 4 * KB, 8 * KB, 16 * KB, 32 * KB, 64 * KB, 128 * KB,
        256 * KB, 512 * KB, 1 * MB, 2 * MB, 4 * MB, 8 * MB
    };
    
    // create and randomize array
    printf_debug("creating array...         ");
    static int arr[ARR_SIZE];
    printf_debug("randomizing array...      ");
    for (int i = 0; i < sizeof(arr)/sizeof(int); i++) {
        arr[i] = (i * 13941) % 35; //13778941
    }
    printf_debug("running tests...         ");
    json results;
    
    
    Timer test_timer;
    test_timer.start();
    cpu_test     (results["cpu"]["reg"], rep_cnt * 100);
    cpu_test_r   (results["cpu"]["read"], arr, sizes, rep_cnt);
    cpu_test_w   (results["cpu"]["write"], arr, sizes, rep_cnt);
    cpu_test_rw  (results["cpu"]["rw"], arr, sizes, rep_cnt);
    io_test_rw   (results["io"]["band"], io_sizes, rep_cnt * 8);
    io_test_many (results["io"]["lat"], (rep_cnt / REP) * KILO);
    test_timer.stop();
    
    printf("---------------------------------\n");
    printf("%-30s: %1.3f\n", "time taken", test_timer.duration.count() * NANO);
    
    printf_debug("generating output...    \n");
    printf("                                ");
    cout << results.dump(true) << endl;
    if (argc >= 2) {
        ofstream ofs (argv[1]);
        ofs << results.dump(true) << endl;
    }
    
    return 0;
}