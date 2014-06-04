/*
 * aplication_base.cc
 *
 */

#ifdef HAVE_PETSC
#include <petsc.h>
#include <petscsys.h>
#endif

#include "system/application_base.hh"
#include "system/sys_profiler.hh"


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

    if ( log_filename == "//" ) {
    	// -l option without given name -> turn logging off
    	sys_info.log=NULL;
    } else	{
    	// construct full log name
    	log_name << log_filename <<  "." << sys_info.my_proc << ".log";
    	sys_info.log_fname = FilePath(log_name.str(), FilePath::output_file );
    	sys_info.log=xfopen(sys_info.log_fname.c_str(),"wt");
    }

    sys_info.verbosity=0;
    sys_info.pause_after_run=0;
}


FILE *ApplicationBase::petsc_output_ =NULL;

#ifdef HAVE_PETSC
PetscErrorCode ApplicationBase::petscvfprintf(FILE *fd, const char format[], va_list Argp) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (fd != stdout && fd != stderr) { /* handle regular files */
    ierr = PetscVFPrintfDefault(fd,format,Argp); CHKERRQ(ierr);
  } else {
    const int buf_size = 65000;
    char buff[65000];
    size_t length;
    ierr = PetscVSNPrintf(buff,buf_size,format,&length,Argp);CHKERRQ(ierr);

    /* now send buff to whatever stream or whatever you want */
    fwrite(buff, sizeof(char), length, petsc_output_);
  }
  PetscFunctionReturn(0);
}
#endif


void ApplicationBase::petsc_initialize(int argc, char ** argv) {
#ifdef HAVE_PETSC
    if (petsc_redirect_file_ != "") {
        petsc_output_ = fopen(petsc_redirect_file_.c_str(), "w");
        PetscVFPrintf = this->petscvfprintf;
    }


    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

    int mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    xprintf(Msg, "MPI size: %d\n", mpi_size);
#endif
}



int ApplicationBase::petcs_finalize() {
#ifdef HAVE_PETSC
	if ( petsc_initialized )
	{
		PetscErrorCode ierr=0;

		ierr = PetscFinalize(); CHKERRQ(ierr);

		if (petsc_output_) fclose(petsc_output_);

		petsc_initialized = false;

		return ierr;
	}
#endif

	return 0;
}


void ApplicationBase::init(int argc, char ** argv) {
    Profiler::initialize();
    // parse our own command line arguments, leave others for PETSc
	this->parse_cmd_line(argc, argv);

	this->petsc_initialize(argc, argv);

    this->system_init(PETSC_COMM_WORLD, log_filename_); // Petsc, open log, read ini file


    this->run();

	this->after_run();
}


ApplicationBase::~ApplicationBase() {
	if (sys_info.log) xfclose(sys_info.log);
	petcs_finalize();
}

