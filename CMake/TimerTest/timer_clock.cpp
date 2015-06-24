#include <stdio.h>
#include <iostream>
#include <time.h>
using namespace std;

const int MEASUREMENTS = 100;

/**
 * Performs minimal time step using 'clock' call
 * @returns value in microseconds
 */
double getResolution () {
    clock_t t1, t2;
    t2 = t1 = clock ();
    while (t1 == t2) {
        t2 = clock ();
    }
    double min_time_step = 1000000.0 * (t2 - t1) / CLOCKS_PER_SEC;
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
