/*
 * aplication_base.hh
 *
 */

#ifndef APPLICATION_BASE_HH_
#define APPLICATION_BASE_HH_


#include <string>
#include <sstream>
#include <mpi.h>

#include "global_defs.h"
#include "system/xio.h"
#include "system/file_path.hh"


using namespace std;


static bool petsc_initialized = false;


/**
 *
 */
class ApplicationBase {
public:

	ApplicationBase(int argc,  char ** argv);

	virtual ~ApplicationBase();

protected:

	virtual void run() = 0;

	virtual void before_flow_123d_init(const int argc, char ** argv) = 0;

	/**
	 * Read system parameters, open log.
	 */
	void system_init( MPI_Comm comm, const string &log_filename );

	/**
	 * Initialize PETSC.
	 */
	void petsc_initialize(int argc, char ** argv);

	/**
	 * Finalize PETSC. If finalization failed return nonzero value.
	 */
	int petcs_finalize();

    /**
     * Log file name argument - passed to system_init; "" menas default, "\n" means no logging
     * TODO: move whole system_init into Application, use singleton for some runtime global options
     * for the Flow123d library.
     */
    string log_filename_;
};

#endif /* APPLICATION_BASE_HH_ */
