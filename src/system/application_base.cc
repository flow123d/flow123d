/*
 * aplication_base.cc
 *
 */

#include <petsc.h>

#include "system/application_base.hh"


ApplicationBase::ApplicationBase(int argc,  char ** argv)
: log_filename_("")
{ }


void ApplicationBase::system_init( MPI_Comm comm, const string &log_filename ) {
    int ierr;

    //for(int i=0;i<argc;i++) xprintf(Msg,"%s,",argv[i]);
    petsc_initialized = true;
    sys_info.comm=comm;


    Xio::init(); //Initialize XIO library

    // TODO : otevrit docasne log file jeste pred ctenim vstupu (kvuli zachyceni chyb), po nacteni dokoncit
    // inicializaci systemu

    ierr=MPI_Comm_rank(comm, &(sys_info.my_proc));
    ierr+=MPI_Comm_size(comm, &(sys_info.n_proc));
    ASSERT( ierr == MPI_SUCCESS,"MPI not initialized.\n");

    // determine logfile name or switch it off
    stringstream log_name;

    if ( log_filename == "\n" ) {
           // -l option without given name -> turn logging off
           sys_info.log=NULL;
    } else
   if (log_filename != "") {
      // given log name
           log_name << log_filename <<  "." << sys_info.my_proc << ".log";
           sys_info.log_fname = FilePath(log_name.str(), FilePath::output_file );
           sys_info.log=xfopen(sys_info.log_fname.c_str(),"wt");

    } else {
        // use default name
        log_name << "flow123."<< sys_info.my_proc << ".log";
        sys_info.log_fname = FilePath(log_name.str(), FilePath::output_file );
        sys_info.log=xfopen(sys_info.log_fname.c_str(),"wt");

    }

    sys_info.verbosity=0;
    sys_info.pause_after_run=0;
}



void ApplicationBase::petsc_initialize(int argc, char ** argv) {
#ifdef HAVE_PETSC
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
#endif
}



int ApplicationBase::petcs_finalize() {
#ifdef HAVE_PETSC
	if ( petsc_initialized )
	{
		PetscErrorCode ierr=0;

		ierr = PetscFinalize(); CHKERRQ(ierr);

		petsc_initialized = false;

		return ierr;
	}
#endif

	return 0;
}


void ApplicationBase::init(int argc, char ** argv) {
    // parse our own command line arguments, leave others for PETSc
	this->parse_cmd_line(argc, argv);

	this->petsc_initialize(argc, argv);

    this->system_init(PETSC_COMM_WORLD, log_filename_); // Petsc, open log, read ini file

	//try {
		this->run();
	//} catch (std::exception & e) {
	//	std::cerr << e.what();
	//	exit( exit_failure );
	//}

	this->after_run();
}


ApplicationBase::~ApplicationBase() {
	if (sys_info.log) xfclose(sys_info.log);
	petcs_finalize();
}

