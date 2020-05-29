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


/// Simple exception for tests that accepts string message.
TYPEDEF_ERR_INFO( EI_Text, std::string);
DECLARE_EXCEPTION(ExcLoggerTest, << EI_Text::val );


// Calls various types of logger messages, allow check logger with different settings.
void logger_messages() {
	unsigned int mesh_size = 150;
	double start_time = 0.5;

	MessageOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements.\n"
			<< "... next line" << "\n" << "... next line" << "\n";
	MessageOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << "\n";

	WarningOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements." << "\n"
			<< "... next line" << "\n";
	WarningOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << "\n";

	LogOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements." << "\n";
	LogOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << "\n";

	DebugOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements." << "\n";
	DebugOut().every_proc() << "Output of process: " << LoggerOptions::get_instance().get_mpi_rank() << "\n";

    // print error message
	try {
        THROW(ExcLoggerTest() << EI_Text("The field has set non-unique names of components.\nPlease set unique names.") );
    } catch (ExcLoggerTest & e) {
    	_LOG( Logger::MsgType::error ) << e.what();
    }

	// necessary for setting next test
	LoggerOptions::get_instance().reset();
}


TEST(Logger, no_log_file) {
	// log file is set to empty
	Profiler::instance();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);
	LoggerOptions::get_instance().set_log_file("");

	logger_messages();
}

TEST(Logger, without_init_log_file) {
	// log file is not set > Log and Debug messages are redirect to screen output
	Profiler::instance();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);

	logger_messages();
}

TEST(Logger, log_file_without_mpi) {
	// MPI is not set > random rank of process is generated
	Profiler::instance();
	LoggerOptions::get_instance().set_log_file("without_mpi");

	logger_messages();
}

TEST(Logger, log_file_with_mpi) {
	// full usage of log
	Profiler::instance();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);
	LoggerOptions::get_instance().set_log_file("with_mpi");

	logger_messages();
}

TEST(Logger, mask_manipulator) {
	Profiler::instance();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);
	LoggerOptions::get_instance().set_log_file("manip");

	MessageOut() << "First message to cout and file.\n"
				 << StreamMask::cout << "Second message only to cout." << std::endl
				 << StreamMask::log << "Third message only to file.\n"
				 << (StreamMask::cout | StreamMask::log) << "Fourth message to cout and file.\n";

	LoggerOptions::get_instance().reset();
}

TEST(Logger, fmt_lib) {
	Profiler::instance();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);
	LoggerOptions::get_instance().set_log_file("");

	int i=1;
	double f=0.2;
	std::string s = "ahoj";
	MessageOut() << fmt::format("int: {} double: {} string: {}\n", i, f, s);
	MessageOut().fmt("int: {} double: {} string: {}\n", i, f, s);

	LoggerOptions::get_instance().reset();
}

TEST(FealAssert, warning) {
	Profiler::instance();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);
	LoggerOptions::get_instance().set_log_file("assert_warn");

	std::string s1 = "feal";
    std::string s2 = "assert";
    FEAL_ASSERT(s1.empty() && s2.empty())(s1)(s2).warning("Strings must be empty.");

    // shorter version of macro - "ASSERT" - is not in conflict with external library
    ASSERT(0).warning("Zero value.");

}
