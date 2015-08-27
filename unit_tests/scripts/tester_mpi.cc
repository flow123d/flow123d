#include <mpi.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <new>          // std::bad_alloc

#include <sys/types.h>
#include <unistd.h>

using namespace std;

/*
 * tester_mpi -t TIME -m MEMORY -p
 * 
 * - run TIME seconds
 * - allocates MEMORY kBytes
 * - use MPI if called with -p
 * - reports MPI size
 */ 
int main(int argc, char * argv[]) {
  int time=0;
  int memory=1;
  bool parallel=false;
  int exit_code=0;
  
  // parse arguments
  for(char ** p=argv+1; p < argv + argc; p++) {
    if ( string(*p) == "-t") {
      p++;
      time=atoi(*p);
    } else
    if ( string(*p) == "-m") {
      p++;
      memory=atoi(*p);
    } else
    if ( string(*p) == "-p" ) 
      parallel=true;
    else 
    if ( string(*p) == "-e" ) {
      p++;
      exit_code=atoi(*p);
    } else {   
      cout << "Unknown parameter: '" << *p << "'." << endl;
    }
  }

  // print PID
  cout << "PID: " << getpid() << endl;
  // test MPI
  int mpi_size=1;
  int mpi_rank=0;
  if (parallel) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  }
  
  if (mpi_rank == 0) {
    cout << "MPI size: " << mpi_size << endl;
    
    // wait time seconds
    int t=0;
    while ( t < time ) {
      cout << "wall time[s]: " << t << endl;
      sleep(1);
      t++;
    }  
    cout << endl;
    
    // allocates    
    char *dummy;
    int total=0, size=memory/2;
    for(; size >0; size=(memory-total)/2) {
        // try catch alloc
        try {
            dummy = new char[size*1024];
        } catch (std::bad_alloc& ba) {
            std::cerr << "bad_alloc caught: " << ba.what() << '\n';
        }

        if (dummy) {
          total+=size;
          cout << "allocated [kB]: " << total << endl;
        } else {
          break;
        }
    }    
  }

  if (parallel)   MPI_Finalize();
  
  exit(exit_code);
}
