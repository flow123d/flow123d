#include <stdio.h>
#include <iostream>
#include <chrono>
using namespace std;

const int MEASUREMENTS = 100;

/**
 * Performs minimal time step using chrono 'high_resolution_clock' method
 * @returns value in microseconds
 * @require c++11
 */
double getResolution () {
    chrono::time_point<chrono::high_resolution_clock> t1, t2;
    t2 = t1 = chrono::high_resolution_clock::now ();
    while (chrono::duration_cast<std::chrono::nanoseconds> (t2 - t1).count () == 0) {
        t2 = chrono::high_resolution_clock::now ();
    }
    double min_time_step = chrono::duration_cast<std::chrono::nanoseconds> (t2 - t1).count ();
//    printf("[%1.3f] ", min_time_step);
    return min_time_step / 1000.0;
}

int main (void) {
    double result = 0;
    for (int i = 0; i < MEASUREMENTS; ++i)
        result += getResolution ();

    printf ("%1.3f", result / MEASUREMENTS);

    return 0;
}
