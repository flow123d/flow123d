// Consumer.cpp : Defines the entry point for the console application.
//

// #include "stdafx.h"
#include <iostream>
#include <string>

#include <chrono>
#include <thread>

#define KB 1024
#define MB (KB*KB)
#define MEMORY_STEPS 10

using namespace std;


/**
 * Function will makes current thread sleep for specific amount of seconds
 * @param long time: wait time in seconds
 */
void limit_time(long time) {
	if (time) {
		cout << "time limit:    " << time << " s" << endl;
		int t = 0;
		while (t < time) {
			cout << "wall time[s]:  " << t << endl;
			this_thread::sleep_for(chrono::seconds(1));
			t++;
		}
		cout << "time limit done" << endl << endl;
	} else {
		this_thread::sleep_for(chrono::hours(1));
	}
}

/**
* Function will allocate specific amount of MB
* @param long memory: amount in MB
*/
void limit_memory(long memory) {
	char *bytes;
	if (memory) {
		cout << "memory limit:  " << memory << " MB" << endl;
		const int ALLOC_SIZE = (memory) / MEMORY_STEPS;
		int m = 0;
		while (m < memory) {
			cout << "allocated[MB]: " << m << endl;
			try {
				bytes = new char[ALLOC_SIZE * MB];
				m += ALLOC_SIZE;
			} catch (const std::bad_alloc& ba) {
				std::cerr << "bad_alloc: " << ba.what() << endl;
				throw ba;
			}
			// work on data to put them into Working set
			for (int j = 0; j < ALLOC_SIZE * MB; j++)
				bytes[j] = 0;
		}
		cout << "allocated[MB]: " << m << endl;
		cout << "memory limit done" << endl << endl;
	}
}


void limit_time_memory(long simul, long time, long memory) {
	const long ALLOC_SIZE = ((double)memory * MB) / simul;
	const long SLEEP_TIME = ((double)time * 1000) / simul;

	cout << "ALLOC_SIZE: " << ALLOC_SIZE << endl;
	cout << "SLEEP_TIME: " << SLEEP_TIME << endl;
	cout << endl;


	cout << "0. step" << endl;
	cout << "allocated[MB]: " << 0 << endl;
	cout << "wall time[s]:  " << 0 << endl;
	cout << endl;

	char * bytes;
	for (int i = 1; i <= simul; i++) {
		// malloc
		try {
			bytes = new char[ALLOC_SIZE];
		} catch (const std::bad_alloc& ba) {
			std::cerr << "bad_alloc: " << ba.what() << endl;
			throw ba;
		}
		// work on data to put them into Working set
		for (int j = 0; j < ALLOC_SIZE; j++)
			bytes[j] = 0;

		// sleep
		this_thread::sleep_for(chrono::milliseconds(SLEEP_TIME));

		cout << i << ". step" << endl;
		cout << "allocated[MB]: " << (i * ALLOC_SIZE) / (float)MB << endl;
		cout << "wall time[s]:  " << (i * SLEEP_TIME) / 1000.0 << endl;
		cout << endl;
	}
}

void usage(string msg) {
	cout << msg << endl;
	cout << "usage: Consumer -t [TIME in s] -m [MEMORY in MB] <-s [NO steps]>" << endl;
	exit(1);
}

int main(int argc, char * argv[]) {
	string tmp;

	long time = 0;
	long memory = 0;
	long simul = 0;
	int time_i = 0;
	int memory_i = 0;

	for (int i = 0; i < argc; ++i) {
		tmp = argv[i];
		if (tmp == "-t" || tmp == "--time") {
			time = stol(string(argv[++i]));
			time_i = i;
		}
		else if (tmp == "-m" || tmp == "--memory") {
			memory = stol(string(argv[++i]));
			memory_i = i;
		}
		else if (tmp == "-s" || tmp == "--simultaneous") {
			simul = stol(string(argv[++i]));
		}
	}

	if (time_i == 0 && memory_i == 0)
		usage("Please specify limits");

	if (!simul) {
		// make limits in order they were declared
		if (time_i < memory_i) {
			limit_time(time);
			limit_memory(memory);
		}
		else {
			limit_memory(memory);
			limit_time(time);
		}
	} else {
		if (time == 0 || memory == 0)
			usage("In simultaneous mode must be both limits greater than 0!");

		limit_time_memory(simul, time, memory);
	}
	cout << "Exiting ..." << endl;

	return 0;
}
