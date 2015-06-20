#include <stdio.h>
#include <iostream>
#include <sys/time.h>
using namespace std;

const int MEASUREMENTS = 100;

/**
 * Performs minimal time step using 'gettimeofday' call
 * @returns value in microseconds
 * @require unix-based platform
 */
double getResolution () {
    timeval t1, t2;
    double min_time_step;
    gettimeofday (&t1, NULL);
    gettimeofday (&t2, NULL);
    while ((t2.tv_usec - t1.tv_usec) == 0) {
        gettimeofday (&t2, NULL);
    }
    min_time_step = (t2.tv_usec - t1.tv_usec);
//    printf("[%1.3f] ", min_time_step);
    return min_time_step;
}

int main (void) {
    double result = 0;
    for (int i = 0; i < MEASUREMENTS; ++i)
        result += getResolution ();

    printf("%1.3f", result / MEASUREMENTS);

    return 0;
}
