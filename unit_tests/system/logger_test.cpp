/*
 * logger_test.cpp
 *
 *  Created on: May 10, 2012
 *      Author: jb
 */


#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "system/logger.hh"
#include "system/sys_profiler.hh"


// Calls various types of logger messages, allow check logger with different settings.
void logger_messages() {
	unsigned int mesh_size = 150;
	double start_time = 0.5;

	MessageOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements.\n"
			<< "... next lines" << "\n" << "... in same message" << std::endl;
	MessageOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << std::endl;

	WarningOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements."
			<< std::endl << "... next line in separate message" << std::endl;
	WarningOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << std::endl;

	LogOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements." << std::endl;
	LogOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << std::endl;

	DebugOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements." << std::endl;
	DebugOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << std::endl;

	// flush screen streams for better but not perfect output
	std::cout << std::flush;
	std::cerr << std::flush;

	// necessary for setting next test
	LoggerOptions::get_instance().reset();
}


TEST(LoggerNoLog, full) {
	Profiler::initialize();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);
	LoggerOptions::get_instance().set_log_file("");

	logger_messages();
}

TEST(LoggerWithoutInitLogFile, full) {
	Profiler::initialize();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);

	logger_messages();
}

TEST(LoggerLogFileWithoutMPI, full) {
	Profiler::initialize();
	LoggerOptions::get_instance().set_log_file("without_mpi");

	logger_messages();
}

TEST(LoggerLogFileWithMPI, full) {
	Profiler::initialize();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);
	LoggerOptions::get_instance().set_log_file("with_mpi");

	logger_messages();
}
