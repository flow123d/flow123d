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
#include "system/logger_options.hh"
#include "system/sys_profiler.hh"


// Calls various types of logger messages, allow check logger with different settings.
void logger_messages() {
	unsigned int mesh_size = 150;
	double start_time = 0.5;

	std::cout << "MessageOut: Start of simulation at time ..." << std::endl;
	MessageOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements.\n"
			<< "... next line" << "\n" << "... next line" << "\n";

	std::cout << "MessageOut: Output of process: ..." << std::endl;
	MessageOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << "\n";


	std::cout << "WarningOut: Start of simulation at time ..." << std::endl;
	WarningOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements." << "\n"
			<< "... next line" << "\n";

	std::cout << "WarningOut: Output of process: ..." << std::endl;
	WarningOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << "\n";

    std::cout << "LogOut: Start of simulation at time ..." << std::endl;
	LogOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements." << "\n";

	std::cout << "LogOut: Output of process: ..." << std::endl;
	LogOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << "\n";

    std::cout << "DebugOut: Start of simulation at time ..." << std::endl;
	DebugOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements." << "\n";

	std::cout << "DebugOut: Output of process: ..." << std::endl;
	DebugOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << "\n";

	std::cout << "=================================================\n" <<std::endl;
	// flush screen streams for better but not perfect output
	//std::cout << std::flush;
	//std::cerr << std::flush;

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

TEST(MaskManipulator, full) {
	Profiler::initialize();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);
	LoggerOptions::get_instance().set_log_file("manip");

	MessageOut() << "First message to cout and file.\n"
				 << StreamMask::cout << "Second message only to cout." << std::endl
				 << StreamMask::log << "Third message only to file.\n"
				 << (StreamMask::cout | StreamMask::log) << "Fourth message to cout and file.\n";

	LoggerOptions::get_instance().reset();
}

TEST(FealAssert, warning) {
	Profiler::initialize();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);
	LoggerOptions::get_instance().set_log_file("assert_warn");

	std::string s1 = "feal";
    std::string s2 = "assert";
    FEAL_ASSERT(s1.empty() && s2.empty())(s1)(s2).warning("Strings must be empty.");

    // shorter version of macro - "ASSERT" - is not in conflict with external library
    ASSERT(0).warning("Zero value.");

}
