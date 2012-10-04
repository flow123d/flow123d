/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @ingroup la
 * @brief	Unified interface to various linear solvers
 * @author	Jan Brezina
 *
 * The only internal (linked) solver is PETSC KPS which is already interface to the number of direct
 * and iterative solvers. Further several external solvers are supported: MATLAB, ISOL (due to Pavel Jiranek)
 */

#include <ctype.h>
#include <strings.h>
#include <petscksp.h>
#include <petscviewer.h>

#include "system/system.hh"
#include "la/distribution.hh"
#include "la/solve.h"
#include "la/linsys.hh"

//#include "profiler.hh"

static void solver_set_type( struct Solver *solver );
static void RunExtern( struct Solver *solver,char *cmdline,void (*write_sys)(struct Solver *), void (*read_sol)(struct Solver *) );
static void clean_directory(void);

// internal solvers
static void solver_petsc(struct Solver *solver);

// drivers for external solvers
// MATLAB
static void write_sys_matlab( struct Solver *solver );
static void write_matlab_linsys(LinSys *mtx,int nz);
static void read_sol_matlab( struct Solver *solver );
// ISOL
static void isol_params_init(ISOL_params *par);
static void write_sys_isol( struct Solver *solver );



//==============================================================================
/*!	@brief	Initialize a solver structure
 *
 * 	Initialize all members form the users options.
 *  Possibly initialize specific solver parameters.
 *
 *  @param[in] solver already allocated structure to be initialized
 */
void solver_init(Solver * solver, Input::AbstractRecord in_rec) {
    double solver_accurancy;

    F_ENTRY;
	if ( solver == NULL ) xprintf(PrgErr,"Structure solver not allocated.\n");

	if (in_rec.type() == Solver::get_input_type_petsc() )  {
	    solver->type = PETSC_SOLVER;
	    solver->params  = Input::Record(in_rec).val<string>("options");
	} else
	if (in_rec.type() == Solver::get_input_type_bddc() ) {
	    solver->type = PETSC_MATIS_SOLVER;
	} else {
	    xprintf(UsrErr,"Unsupported solver: %s\n", in_rec.type().type_name().c_str());
	}

	solver->r_tol = Input::Record(in_rec).val<double>("r_tol");
	solver->a_tol = Input::Record(in_rec).val<double>("a_tol");

/*   	solver->executable = OptGetStr( "Solver", "Solver_executable",solver->name);

    solver->keep_files    = OptGetBool( "Solver", "Keep_solver_files", "no" );
    solver->manual_run     = OptGetBool( "Solver", "Manual_solver_run", "no" );
    solver->use_ctrl_file  = OptGetBool( "Solver", "Use_control_file", "no" );
    if (solver->use_ctrl_file)
    	solver->ctrl_file      = IONameHandler::get_instance()->get_input_file_name(OptGetStr( "Solver", "Control_file", NULL )).c_str();
    solver->use_last_sol    =OptGetBool( "Solver", "Use_last_solution", "no" );
    /// Last solution reuse is possible only for external solvers
    if (solver->use_last_sol && (solver->type == PETSC_SOLVER)) {
        xprintf(UsrErr,"Can not reuse last solution, when using an internal solver.");
    }

    //! generic solver parameters
    solver_accurancy=   OptGetDbl("Solver","Solver_accurancy","1.0e-7");
    solver->max_it=     OptGetInt("Solver", "max_it", "200" );
    solver->r_tol=      OptGetDbl("Solver", "r_tol", "-1" );
    if (solver->r_tol < 0) solver->r_tol=solver_accurancy;
    solver->a_tol=      OptGetDbl("Solver", "a_tol", "1.0e-9" );

	if (solver->type == ISOL) {
    	solver->isol_params=(ISOL_params *)malloc(sizeof(ISOL_params));
    	isol_params_init(solver->isol_params);
    }*/
}



Input::Type::AbstractRecord & Solver::get_input_type() {
    using namespace Input::Type;
    static AbstractRecord rec("Solver", "Solver setting.");

    if (!rec.is_finished()) {
        rec.declare_key("a_tol", Double(0.0), Default("1.0e-9"),
                "Absolute residual tolerance.");
        rec.declare_key("r_tol", Double(0.0, 1.0), Default("1.0e-7"),
                "Relative residual tolerance (to initial error).");
        rec.declare_key("max_it", Integer(0), Default("10000"),
                "Maximum number of outer iterations of the linear solver.");
        rec.finish();

        Solver::get_input_type_petsc();
        Solver::get_input_type_bddc();

        rec.no_more_descendants();
    }
    return rec;
}

Input::Type::Record & Solver::get_input_type_petsc() {
    using namespace Input::Type;
    static Record rec("Petsc", "Solver setting.");

    if (!rec.is_finished()) {
        rec.derive_from(Solver::get_input_type());
        rec.declare_key("options", String(), Default(""),  "Options passed to the petsc instead of default setting.");
        rec.finish();
    }
    return rec;
}

Input::Type::Record & Solver::get_input_type_bddc() {
    using namespace Input::Type;
    static Record rec("Bddc", "Solver setting.");

    if (!rec.is_finished()) {
        rec.derive_from(Solver::get_input_type());
        rec.finish();
    }
    return rec;
}


//=============================================================================
/*! @brief  Set solver type from its name.
 *
 *  The name comparison is case insensitive.
 */
/// one test macro
#define TEST_TYPE(name_str,id) if ( strcmpi( solver->name, name_str ) == 0 ) solver->type = (id);

void solver_set_type( Solver *solver )
{
    F_ENTRY;
    solver->type=UNKNOWN;
    TEST_TYPE("petsc",PETSC_SOLVER);
    TEST_TYPE("petsc_matis",PETSC_MATIS_SOLVER);
    TEST_TYPE("si2",SI2);
    TEST_TYPE("gi8",GI8);
    TEST_TYPE("isol",ISOL);
    TEST_TYPE("matlab",MATLAB);

}

//=============================================================================
/*! @brief Solves a given linear system.
 *
 * Call user selected internal or external solver.
 * @param[in] solver solver structure to use
 * @param[in,out] system linear system to solve, conatains also result
 *
 * @pre Initialized solver. Assembled system.
 * @post Valid solution.
 */

void solve_system( struct Solver *solver, struct LinSys *system )
{
/// set command line for external solvers
#define SET_GENERIC_CALL sprintf( cmdline, "%s %s",solver->executable,solver->params.c_str())
#define SET_MATLAB_CALL sprintf( cmdline, "matlab -r solve" )

	char cmdline[ 1024 ];

	ASSERT( NONULL( solver ),"NULL 'solver' argument.\n");
    ASSERT( NONULL( system ),"NULL 'system' argument.\n");

	solver->LinSys=system;
	switch (solver->type) {
			// internal solvers
	        case PETSC_SOLVER:
	        case PETSC_MATIS_SOLVER:
	        	solver_petsc( solver );
	            break;
	        // external solvers
	      	case ISOL:
	      		SET_GENERIC_CALL;
	      		RunExtern(solver,cmdline,&write_sys_isol, &read_sol_matlab);
	      	case MATLAB:
	      		SET_MATLAB_CALL;
	      		RunExtern(solver,cmdline,&write_sys_matlab,&read_sol_matlab);
	        	break;
	      	case GI8:
	        case SI2:
	            xprintf(UsrErr,"Solver %s is not supported.\n",solver->name);
	        case UNKNOWN:
	            xprintf(UsrErr,"UNKNOWN solver is not supported.\n");
	}

}

//=============================================================================
/*! @brief Call an external solver.
 *
 *  Make temporary directory, write down matrix and RHS vector. Then
 *  perform system call of given program and read solution form given file.
 *
 *  @param[in,out] solver solver structure to use (solution in linear system)
 *  @param[in] cmdline calling command line
 *  @param[in] write_sys function to write down linear system
 *  @param[in] read_sol function to read the solution
 */


void RunExtern( Solver *solver,char *cmdline,
		void (*write_sys)(Solver *), void (*read_sol)(Solver *) )
{
	char cmd[ LINE_SIZE ];

	xmkdir( "tmp" );
	xchdir( "tmp" );
	// TODO: use_last_sol by melo mit parametr s odkazem kde jsou data
	// a jakeho typu, malo by jit ulozit i z vnitrniho solveru
	if (! solver->use_last_sol) {
	    /// prepare data
	    clean_directory();
	    xprintf( Msg, "Writing files for solver... ")/*orig verb 2*/;
	    (*write_sys)(solver); // prepare solver specific input files
	    xprintf(Msg, "run extr. solver\n");
	    /// if user provides an explicit control file, use that one instead
	    if( solver->use_ctrl_file == true ) {
			sprintf( cmd, "copy /Y ../%s ", solver->ctrl_file );
			xsystem( cmd );
	    }

	    /// run the solver
	    if( solver->manual_run == true ) {
	        xprintf( Msg, "\nPROGRAM PAUSED.\n"
				"Start solver of linear equations manually:\n\n%s\n\n"
				"Press ENTER when calculation finished.\n" );
	        getchar();
	    } else {
	        xprintf( Msg, "\nCalling solver... \n"
					 "BEGIN OF SOLVER'S MESSAGES\n");
	        strcpy(cmd,cmdline);
	        // pipe stdout and stderr through tee to get it to both, the screen and the logfile
	        strcat(cmd," 2>&1 | tee ../flow_extern_solver.log");
	        xsystem( cmd );
	    }
	}
	/// read the solution
	(*read_sol)( solver );

	if( solver->keep_files == false ) {
		clean_directory();
		xchdir( ".." );
		xremove( "tmp" );
	} else {
		xchdir( ".." );
	}

    xprintf( Msg, "END OF MESSAGES OF THE SOLVER\n\n");
}
//=============================================================================
/*! @brief Clean temporary directory of external solver.
 */

void clean_directory( void )
{
	xremove( "rid.ghp" );
	xremove( "vector.sro" );
	xremove( "vector.sri" );
	xremove( "vector.log" );
	xremove( "solve.m" );
	xremove( "matrix.dat" );
	xremove( "rhs.dat" );
	xremove( "solution.dat" );
}

//=========================================================================
/*! @brief  solve given system by internal PETSC KSP solver.
 *
 * - insert given or default option string into the PETSc options db
 * - create KSP
 * - set tolerances
 * - solve, report number of iterations and convergency reason
 *
 * default choice of PETSC solver options:
 *
 *  - for SPD system (using NSchurs=1,2):\n
 *          "-ksp_type cg -pc_type ilu -pc_factor_levels 5 -ksp_diagonal_scale_fix"\n
 *          (CG solver and ILU preconditioner with 5 levels and diagonal scaling)
 *  - for other system (NSchurs=0):\n
 *          "-ksp_type bcgs -pc_type ilu -pc_factor_levels 5 -ksp_diagonal_scale fix"\n
 *          (BiCGStab solver and ILU preconditioner with 5 levels and diagonal scaling)\n
 *
 * For bad systems you can try:
 *          - remove diagonal scaling
 *          - use additional option: "-pc_factor_shift_nonzero"
 *
 * To this end you can use either the command line or <b> [Solver] petsc_options </b> parameter.
 *
 *
 * @}
 */

void solver_petsc(Solver *solver)
{
	LinSys *sys=solver->LinSys;
	KSP System;
	KSPConvergedReason Reason;
        //PetscViewer mat_view;
	const char *petsc_dflt_opt;
	const char *petsc_str;
	int nits;

	F_ENTRY;

	//LSView(sys);

	if (solver->type == PETSC_MATIS_SOLVER) {
           if (sys->ds().np() > 1) {

	       // parallel setting
              if (sys->is_positive_definite())
                  petsc_dflt_opt="-ksp_type cg -pc_type nn -nn_coarse_pc_factor_mat_solver_package mumps -is_localD_pc_factor_mat_solver_package mumps -is_localN_pc_factor_mat_solver_package mumps";
              else
                  if (sys->is_symmetric())
                     petsc_dflt_opt="-ksp_type minres -pc_type nn -nn_coarse_pc_factor_mat_solver_package mumps -is_localD_pc_factor_mat_solver_package mumps -is_localN_pc_factor_mat_solver_package mumps";
	          else
                     petsc_dflt_opt="-ksp_type bcgs -pc_type nn -nn_coarse_pc_factor_mat_solver_package mumps -is_localD_pc_factor_mat_solver_package mumps -is_localN_pc_factor_mat_solver_package mumps";

	   } else {
	       // serial setting
              if (sys->is_positive_definite())
                  petsc_dflt_opt="-ksp_type cg -pc_type nn -nn_coarse_pc_factor_mat_solver_package mumps -is_localD_pc_factor_mat_solver_package mumps -is_localN_pc_factor_mat_solver_package mumps";
              else
                  if (sys->is_symmetric())
                     petsc_dflt_opt="-ksp_type minres -pc_type nn -nn_coarse_pc_factor_mat_solver_package mumps -is_localD_pc_factor_mat_solver_package mumps -is_localN_pc_factor_mat_solver_package mumps";
	          else
                     petsc_dflt_opt="-ksp_type bcgs -pc_type nn -nn_coarse_pc_factor_mat_solver_package mumps -is_localD_pc_factor_mat_solver_package mumps -is_localN_pc_factor_mat_solver_package mumps";
	   }
	}
	else
	{
	   // -mat_no_inode ... inodes are usefull only for
           //  vector problems e.g. MH without Schur complement reduction	
           if (sys->ds().np() > 1) {
	       // parallel setting
              if (sys->is_positive_definite())
                  petsc_dflt_opt="-ksp_type cg -ksp_diagonal_scale_fix -pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 3 -sub_pc_factor_shift_positive_definite -sub_pc_factor_fill 6.0";
                  //petsc_dflt_opt="-ksp_type preonly -pc_type cholesky -pc_factor_mat_solver_package mumps -mat_mumps_sym 1";
                  // -ksp_type preonly -pc_type lu 
              else
                  petsc_dflt_opt="-ksp_type bcgs -ksp_diagonal_scale_fix -pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 3";

	   } else {
	       // serial setting
              if (sys->is_positive_definite())
                  petsc_dflt_opt="-ksp_type cg -pc_type ilu -pc_factor_levels 3 -ksp_diagonal_scale_fix -pc_factor_shift_positive_definite -pc_factor_fill 6.0";
              else
                  petsc_dflt_opt="-ksp_type bcgs -pc_type ilu -pc_factor_levels 5 -ksp_diagonal_scale_fix";
	   }
	}

	if (solver->params == "") petsc_str=petsc_dflt_opt;
	else petsc_str=solver->params.c_str();
	        //OptGetStr("Solver","Solver_params",petsc_dflt_opt);

	xprintf(MsgVerb,"inserting petsc options: %s\n", petsc_str );
	 
	PetscOptionsInsertString(petsc_str); // overwrites previous options values
	//xfree(petsc_str);
    
        MatSetOption(sys->get_matrix(), MAT_USE_INODES, PETSC_FALSE);


    //    xprintf(Msg,"View KSP system\n");
        //Mat matrixForPrint;
        //PetscErrorCode ierr;
        //PetscInt m, n;
        //MatGetSize( sys->get_matrix(), &m, &n );

        //ierr = MatCreate( PETSC_COMM_WORLD, &matrixForPrint ); CHKERRV( ierr ); 
        //ierr = MatSetType( matrixForPrint, MATMPIAIJ ); CHKERRV( ierr ); 
        //ierr = MatSetSizes( matrixForPrint, PETSC_DECIDE, PETSC_DECIDE, m, n ); 

        //std::cout << "Size of the matrix is :" << m << ", " << n << std::endl;
        //Vec auxIn, auxOut;
        //for ( int i = 0; i < n; i++ ) {
        //    // create auxiliary vector of unit matrix
        //    ierr = MatGetVecs( sys->get_matrix(), &auxIn, &auxOut ); CHKERRV( ierr );

        //    VecSetValue( auxIn, i, 1., INSERT_VALUES );
        //    ierr = VecAssemblyBegin( auxIn ); CHKERRV( ierr ); 
        //    ierr = VecAssemblyEnd(   auxIn ); CHKERRV( ierr ); 
 
        //    ierr = MatMult( sys->get_matrix(), auxIn, auxOut ); CHKERRV( ierr ); 

        //    PetscInt low, high;
        //    VecGetOwnershipRange( auxOut, &low, &high );
        //    PetscInt locSize = high - low;

        //    PetscScalar *values;
        //    VecGetArray( auxOut, &values );

        //    std::vector<PetscInt> rows;
        //    std::vector<PetscInt> columns;
        //    for ( int j = low; j < high; j++ ) {
        //        rows.push_back( j );
        //    }
        //    columns.push_back( i );

        //    MatSetValues( matrixForPrint, locSize, &(rows[0]), 1, &(columns[0]), values, INSERT_VALUES );

        //    VecRestoreArray( auxOut, &values );
        //    VecDestroy( auxIn );
        //    VecDestroy( auxOut );
        //}
        //ierr = MatAssemblyBegin( matrixForPrint, MAT_FINAL_ASSEMBLY ); CHKERRV( ierr ); 
        //ierr = MatAssemblyEnd(   matrixForPrint, MAT_FINAL_ASSEMBLY ); CHKERRV( ierr ); 


        //PetscViewer matViewer;
        //PetscViewerASCIIOpen( PETSC_COMM_WORLD, "matrix.m", &matViewer );
        //PetscViewerSetFormat(matViewer,PETSC_VIEWER_ASCII_MATLAB);
        //MatView( matrixForPrint, matViewer );
        //MatDestroy( matrixForPrint );
        //PetscViewerDestroy(matViewer);

        //PetscViewer rhsViewer;
        //PetscViewerASCIIOpen( PETSC_COMM_WORLD, "rhs.m", &rhsViewer );
        //PetscViewerSetFormat(rhsViewer,PETSC_VIEWER_ASCII_MATLAB);
        //VecView( sys->get_rhs(), rhsViewer );
        //PetscViewerDestroy(rhsViewer);

	KSPCreate(PETSC_COMM_WORLD,&System);
	KSPSetOperators(System, sys->get_matrix(), sys->get_matrix(), DIFFERENT_NONZERO_PATTERN);
	KSPSetTolerances(System, solver->r_tol, solver->a_tol, PETSC_DEFAULT,PETSC_DEFAULT);
	KSPSetFromOptions(System);
	KSPSolve(System, sys->get_rhs(), sys->get_solution());
	KSPGetConvergedReason(System,&Reason);
	KSPGetIterationNumber(System,&nits);

	// TODO: make solver part of LinSyt, and make gatter for num of it
	xprintf(MsgLog,"convergence reason %d, number of iterations is %d\n", Reason, nits);
    Profiler::instance()->set_timer_subframes("SOLVING MH SYSTEM", nits);
	KSPDestroy(&System);

}

//=============================================================================
//=============================================================================
/// @brief ISOL calling functions
//!@{

//=============================================================================
/* @brief Read ISOL specific parameters.
 *
 */
void isol_params_init(ISOL_params *par) {

        F_ENTRY;
        /*
        par->method           = OptGetStr( "Solver parameters", "method", "fgmres" );
        par->restart          = OptGetInt( "Solver parameters", "restart", "20");
        par->stop_crit        = OptGetStr( "Solver parameters", "stop_crit", "backerr" );
        par->be_tol           = OptGetDbl( "Solver parameters", "be_tol", "1e-10" );
        par->stop_check       = OptGetInt( "Solver parameters", "stop_check", "1" );
        par->scaling          = OptGetStr( "Solver parameters", "scaling", "mc29_30" );
        par->precond          = OptGetStr( "Solver parameters", "precond", "ilu" );
        par->sor_omega        = OptGetDbl( "Solver parameters", "sor_omega", "1.0" );
        par->ilu_cpiv         = OptGetInt( "Solver parameters", "ilu_cpiv", "0" );
        par->ilu_droptol      = OptGetDbl( "Solver parameters", "ilu_droptol", "1e-3" );
        par->ilu_dskip        = OptGetInt( "Solver parameters", "ilu_dskip", "-1" );
        par->ilu_lfil         = OptGetInt( "Solver parameters", "ilu_lfil", "-1" );
        par->ilu_milu         = OptGetInt( "Solver parameters", "ilu_milu", "0" );
        */
}

//=============================================================================
/*! @brief  Write down input for ISOL
 *
 */
void write_sys_isol( struct Solver *solver )
{
	FILE *         fconf;
    ISOL_params *   par=solver->isol_params;

    F_ENTRY;
    /// write the control file
    fconf = xfopen( "isol.conf", "wt" );
    xfprintf( fconf, "mat_file  = matrix.dat\n" );
    xfprintf( fconf, "mat_fmt   = coo\n" );
    xfprintf( fconf, "rhs_file  = rhs.dat\n" );
    xfprintf( fconf, "mat_sym   = 1\n" );
    xfprintf( fconf, "mat_id_off    = 0\n" );
    xfprintf( fconf, "mat_symmetrize  = 1\n" );
    xfprintf( fconf, "sol0_file  = 0\n" );
    xfprintf( fconf, "sol_file   = solution.dat\n" );
    DBGMSG("first\n");
    xfprintf( fconf, "method = %s\n", par->method);
    DBGMSG("second\n");
    xfprintf( fconf, "max_it   = %d\n",solver->max_it );
    xfprintf( fconf, "max_dim  = %d\n",par->restart );
    xfprintf( fconf, "stop_crit = %s\n",par->stop_crit );
    xfprintf( fconf, "rel_tol = %e\n",solver->r_tol );
    xfprintf( fconf, "abs_tol = %e\n", solver->a_tol);
    xfprintf( fconf, "be_tol  = %e\n",par->be_tol );
    xfprintf( fconf, "stop_check = %d\n",par->stop_check );
    xfprintf( fconf, "scaling = %s\n",par->scaling );
    xfprintf( fconf, "precond = %s\n",par->precond );
    xfprintf( fconf, "sor_omega = %e\n",par->sor_omega );
    xfprintf( fconf, "ilu_droptol = %e\n",par->ilu_droptol );
    xfprintf( fconf, "ilu_milu    = %d\n",par->ilu_milu );
    xfprintf( fconf, "ilu_cpiv    = %d ! do not use at the moment, needs some fixing!\n",
            par->ilu_cpiv );
    if(par->ilu_dskip == -1)
              xprintf(PrgErr,"Sorry, ilu_dskip parameter is not supported.\n");
//            par->ilu_dskip = solver->LinSys->sizeA;
    xfprintf( fconf, "ilu_dskip   = %d\n",par->ilu_dskip );
    xfprintf( fconf, "ilu_lfil    = %d\n",par->ilu_lfil);
    xfclose( fconf );

    /// write the linear system - using MATLAB format
    write_matlab_linsys(solver->LinSys,1);
}
//!@}
// ISOL block

//void write_sys_gm6( Solver *solver )
//void write_vector_sro( LinSystem *mtx )
//void read_sol_gm6( Solver *solver )
//============================================================================
//============================================================================
/// @brief MATLAB call functions
//!@{

//=============================================================================
/*! @brief Write input for MATLAB
 *
 *  written files:  matrix.dat, rhs.dat, solve.m
 */
void write_sys_matlab( struct Solver *solver )
{
	FILE *out;

	//solve.m
	out = xfopen( "solve.m", "wt" );
	xfprintf( out, "load matrix.dat\n" );
	xfprintf( out, "a = spconvert( matrix )\n" );
//	xfprintf( out, "a = a + tril( transpose( a ), -1 )\n" ); // be sure about symmetry
	xfprintf( out, "load rhs.dat\n" );
	xfprintf( out, "sol = a \\ rhs\n" );
	xfprintf( out, "save solution.dat sol -ascii -double\n" );
	xfprintf( out, "exit\n" );
	xfclose( out );

	write_matlab_linsys(solver->LinSys,0);
}

//==========================================================================
/*! @brief Write MATLAB system
 *
 *  Write matrix in MATLAB format to matrix.dat and RHS into rhs.dat.
 *  @param[in] mtx Linear system to write out.
 *  @param[in] write_nz (1 - write number of non-zeroes (ISOL); 0 - don't write (MATLAB))
 */
void write_matlab_linsys(LinSys *mtx,int write_nz) {
    //FILE *out;
    //int mi, ji,nnz;

    F_ENTRY;

    // use LinSys view to write down system in MATLAB format
/*
    LSSetCSR( mtx );
    if (write_nz) nnz=mtx->i[mtx->size];
    else nnz=0;
    /// matrix.dat
    out = xfopen( "matrix.dat", "wt" );
    xfprintf( out, "%d %d %d\n", mtx->size, mtx->size, nnz);
    ji = 0;
    for( mi = 1; mi <= mtx->size; mi++ )
        for( ; ji < mtx->i[mi] ; ji++ )
            xfprintf( out, "%d %d %.16lg\n",mi,mtx->j[ji] + 1, mtx->a[ ji ] );
    xfclose( out );
    /// rhs.dat
    out = xfopen( "rhs.dat", "wt" );
    for( mi = 0; mi < mtx->size; mi++ ) xfprintf( out, "%.16lg\n", mtx->vb[ mi ] );
    xfclose( out );
    */
}

//=============================================================================
/*! @brief Read MATLAB solution into solver->LinSys->vx
 */
void read_sol_matlab( struct Solver *solver )
{
	LinSys *sys=solver->LinSys;
	FILE *in;
	double value;
	int mi;

	in = xfopen( "solution.dat", "rt" );
	int loc_row=0;
	for( mi = 0; mi < sys->size(); mi++ )
	{
        xfscanf( in, "%lf", value );
        if (sys->ds().is_local(mi)) *(sys->get_solution_array() + loc_row)=value;
	}
	xfclose( in );
}

//-----------------------------------------------------------------------------
// vim: set cindent:
