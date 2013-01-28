/*
 * backtrace_test.cpp
 *
 *  Created on: Jan 28, 2013
 *      Author: jb
 */



#include <gtest/gtest.h>
#include <execinfo.h>
#include <iostream>
#include <cstdlib>

using namespace std;

TEST(Backtrace, all) {
    const unsigned int size=1000;
    void ** buffer = new void * [size];
    char  **output;

    int n_entries= backtrace(buffer, size);
    output = backtrace_symbols (buffer, n_entries);
    cout << "Backtrace:" << endl;
    for(unsigned int i=0; i<n_entries; i++)
        cout << output[i] << endl;
    cout << "end of Backtrace" << endl;

    delete [] buffer;
    free(output);


}
