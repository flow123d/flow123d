#define TEST_USE_PETSC

#include "flow_gtest_mpi.hh"
#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "la/distribution.hh"
#include "la/local_system.hh"

#include <petscmat.h>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include <armadillo>
#include "mpi.h"


const int 
    m = 9,
    n_cols = 2,
    n_loops = 5e6,
    ls_size = 100,
    offset = 20;


void allocate_linsys(LinSys* ls){
    ls->start_allocation();
    for(unsigned int i = 0; i<ls_size; i++)
        for(unsigned int j = 0; j<ls_size; j++)
            ls->mat_set_value(i,j,0.0);
}

void print_matrix(std::string matlab_file, const Mat* mat){
    std::string output_file = FilePath(matlab_file + ".m", FilePath::output_file);
    PetscViewer    viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, output_file.c_str(), &viewer);
    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    MatView( *const_cast<Mat*>(mat), viewer);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(LinSys_PETSC, mat_set_value) {
    // setup FilePath directories
    FilePath::set_io_dirs(".",string(UNIT_TESTS_SRC_DIR)+"/la/","",".");
    
    LinSys * ls = new LinSys_PETSC(new Distribution(ls_size, MPI_COMM_WORLD));
    allocate_linsys(ls);
    ls->start_add_assembly();
    
    int i,j;
    int from = offset,
        to = offset+m;
    START_TIMER("LinSys_PETSC_mat_set_value_mm");
    
    for(unsigned int q = 0; q < n_loops; q++)
        for(i = from; i<to; i++)
            for(j = from; j<to; j++)
                ls->mat_set_value(i, j, 1.0);
    
    END_TIMER("LinSys_PETSC_mat_set_value_mm");
    EXPECT_TIMER_LE("LinSys_PETSC_mat_set_value_mm", 22);
    
    ls->finish_assembly();
//     print_matrix("ls_matrix", ls->get_matrix());
    
//     static ofstream os( FilePath("pr_mat_set_value_mm.log", FilePath::output_file) );
//     Profiler::instance()->output(MPI_COMM_WORLD, os);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(LinSys_PETSC, mat_set_values_mm) {
    // setup FilePath directories
    FilePath::set_io_dirs(".",string(UNIT_TESTS_SRC_DIR)+"/la/","",".");
    
    LinSys * ls = new LinSys_PETSC(new Distribution(ls_size, MPI_COMM_WORLD));
    allocate_linsys(ls);
    ls->start_add_assembly();
    
    // 'local system' to be added
    double vals[m*m];
    for(int i = 0; i<m*m; i++) vals[i] = 1.0;
    int rows[m], cols[m];
    for(int i = 0; i<m; i++){
        rows[i] = offset + i;
        cols[i] = offset + i;
    }
        
    
    START_TIMER("LinSys_PETSC_mat_set_values_mm");
    
    for(unsigned int q = 0; q < n_loops; q++)
        ls->mat_set_values(m, rows, m, cols, vals);
    
    END_TIMER("LinSys_PETSC_mat_set_values_mm");
    EXPECT_TIMER_LE("LinSys_PETSC_mat_set_values_mm", 4.5);
    
    ls->finish_assembly();
//     print_matrix("ls_matrix", ls->get_matrix());
    
//     static ofstream os( FilePath("pr_mat_set_values_mm.log", FilePath::output_file) );
//     Profiler::instance()->output(MPI_COMM_WORLD, os);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(LinSys_PETSC, mat_set_values_mcols) {
    // setup FilePath directories
    FilePath::set_io_dirs(".",string(UNIT_TESTS_SRC_DIR)+"/la/","",".");
    
    LinSys * ls = new LinSys_PETSC(new Distribution(ls_size, MPI_COMM_WORLD));
    allocate_linsys(ls);
    ls->start_add_assembly();
    
    // 'local system' to be added
    double vals[m*n_cols];
    for(int i = 0; i<m*n_cols; i++) vals[i] = 1.0;
    int rows[m], cols[n_cols];
    for(int i = 0; i<m; i++) rows[i] = offset + i;
    for(int i = 0; i<n_cols; i++) cols[i] = offset + i;
    
        
    
    START_TIMER("LinSys_PETSC_mat_set_values_mcols");
    
    for(unsigned int q = 0; q < n_loops; q++)
        ls->mat_set_values(m, rows, n_cols, cols, vals);
    
    END_TIMER("LinSys_PETSC_mat_set_values_mcols");
    EXPECT_TIMER_LE("LinSys_PETSC_mat_set_values_mcols", 1.5);
    
    ls->finish_assembly();
//     print_matrix("ls_matrix", ls->get_matrix());
    
//     static ofstream os( FilePath("pr_mat_set_values_mcols.log", FilePath::output_file) );
//     Profiler::instance()->output(MPI_COMM_WORLD, os);
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(LinSys_PETSC, set_local_system_mm) {
    // setup FilePath directories
    FilePath::set_io_dirs(".",string(UNIT_TESTS_SRC_DIR)+"/la/","",".");
    
    LinSys * ls = new LinSys_PETSC(new Distribution(ls_size, MPI_COMM_WORLD));
    allocate_linsys(ls);
    ls->start_add_assembly();
    
    // 'local system' to be added
    LocalSystem loc(m,m);
    for(int i = 0; i<m; i++){
        loc.row_dofs[i] = offset + i;
        loc.col_dofs[i] = offset + i;
        for(int j = 0; j<m; j++)
            loc.add_value(i,j,1.0);
    }
        
    START_TIMER("LinSys_PETSC_set_local_system_mm");
    
    for(unsigned int q = 0; q < n_loops; q++)
        ls->set_local_system(loc);
    
    END_TIMER("LinSys_PETSC_set_local_system_mm");
    EXPECT_TIMER_LE("LinSys_PETSC_set_local_system_mm", 6);
    
    ls->finish_assembly();
//     print_matrix("ls_matrix", ls->get_matrix());
    
//     static ofstream os( FilePath("pr_set_local_system_mm.log", FilePath::output_file) );
//     Profiler::instance()->output(MPI_COMM_WORLD, os);
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(PETSC_mat, mat_set_values_mm) {
    // setup FilePath directories
    FilePath::set_io_dirs(".",string(UNIT_TESTS_SRC_DIR)+"/la/","",".");
    
    PetscErrorCode ierr;
    Mat mat;
    
    ierr = MatCreate(MPI_COMM_WORLD, &mat); CHKERRV( ierr );
    ierr = MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, ls_size, ls_size); CHKERRV( ierr );
    ierr = MatSetType(mat,MATAIJ); CHKERRV( ierr );
    ierr = MatSetUp(mat); CHKERRV( ierr );
    ierr = MatAssemblyBegin(mat, MAT_FLUSH_ASSEMBLY); CHKERRV( ierr );
    
    // 'local system' to be added
    double vals[m*m];
    for(int i = 0; i<m*m; i++) vals[i] = 1.0;
    int rows[m], cols[m];
    for(int i = 0; i<m; i++){
        rows[i] = offset + i;
        cols[i] = offset + i;
    }
    
    
    START_TIMER("PETSC_mat_petsc_mat_set_values");
    
    for(unsigned int q = 0; q < n_loops; q++)
        ierr = MatSetValues(mat, m, rows, m, cols, vals, ADD_VALUES); CHKERRV( ierr );
    
    END_TIMER("PETSC_mat_petsc_mat_set_values");
    EXPECT_TIMER_LE("PETSC_mat_petsc_mat_set_values", 4.5);
    
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY); CHKERRV( ierr );
//     print_matrix("petsc_matrix", &mat);
    
//     static ofstream os( FilePath("pr_petsc_mat_set_values.log", FilePath::output_file) );
//     Profiler::instance()->output(MPI_COMM_WORLD, os);
    static ofstream os( FilePath("set_values_profiler.log", FilePath::output_file) );
    Profiler::instance()->output(MPI_COMM_WORLD, os);
}
