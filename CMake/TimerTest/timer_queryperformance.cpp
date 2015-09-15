#include <stdio.h>
#include <iostream>
#include <windows.h>
using namespace std;

const int MEASUREMENTS = 100;

/**
 * Performs minimal time step using 'QueryPerformanceFrequency' call
 * @returns value in microseconds
 * @require windows platform
 */
double getResolution () {
    LARGE_INTEGER frequency;        // ticks per second
    LARGE_INTEGER t1, t2;           // ticks
    double min_time_step;
    QueryPerformanceFrequency (&frequency);
    QueryPerformanceCounter (&t1);
    QueryPerformanceCounter (&t2);
    min_time_step = (t2.QuadPart - t1.QuadPart) * 1000000.0 / frequency.QuadPart;
    while (min_time_step == 0) {
        QueryPerformanceCounter (&t2);
        min_time_step = (t2.QuadPart - t1.QuadPart) * 1000000.0 / frequency.QuadPart;
    }
//    printf("[%1.3f] ", min_time_step);
    return min_time_step;
}

int main (void) {
    double result = 0;
    for (int i = 0; i < MEASUREMENTS; i++)
        result += getResolution ();

    printf("%1.3f", result / MEASUREMENTS);

    return 0;
}
