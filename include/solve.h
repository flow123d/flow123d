/*!
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @brief Unified interface to various solvers
 *
 * aim: simplified interface to PETSC solvers, Mat, Vec
 * - independent of Seq-MPI
 * - application to parallel PETSC mapping
 * - preallocation by abstract assembly
 *
 * TODO:
 * - samostatny objekt LS
 * - samostatne objekty pro matici a vektor
 * - objekty pro mapovani indexu a jejich automaticke pouziti pri asemblaci
 * - podminene veci jako LS View a old_4_new
 * - funkcni stavy LS
 * - solve jako metoda
 * - nastaveni parametru na volajici strane (zavisle na problemu)
 */

#ifndef SOLVE_H
#define SOLVE_H

/*!************************************************
 * ISOL specific parameters
 **************************************/
// Specific parameters for ISOL solver
typedef struct ISOL_params {
    char*           method;         //!< Iteration method
    int             restart;        //!< num of iter. of restart of GMRES
    char*           stop_crit;      //!< Stoping criterion
    double          be_tol;         //!< Backward error tolerance
    int             stop_check;     //!< Stop check
    char*           scaling;        //!< Scaling method
    char*           precond;        //!< Type od preconditioning
    double          sor_omega;      //!< Relaxation parameter
    /// @name ILU preconditioner parameters
    ///@{
    double          ilu_droptol;      //!< drop tolerance
    int             ilu_milu;         //
    int             ilu_cpiv;
    int             ilu_dskip;
    int             ilu_lfil;
    ///@}
} ISOL_params;


/*!************************************************************
 * @name Solver Type
 */
typedef enum {
    UNKNOWN=0,
    SI2=1,
    GI8=2,
    MATLAB=4,
    PETSC_SOLVER=6,
    ISOL=7,
    PETSC_MATIS_SOLVER=8
} SolverType;


/*!**************************************************************
 *  Solver structure - no matrix, but parameters ...
 **************************************************************/
typedef struct Solver {
	SolverType		type;	  		//!< type of the solver
	char	*name;   		//!< Name of the solver
	char	*executable;  	//!< full path to the external solver executable file
	struct 	LinSys   *LinSys;  	//!< System to solve
	char    external;      	//!< run an external progam as a solver
	char	manual_run;		//!< Run solver manualy ?
	char	use_ctrl_file;	//!< User provided control file ?
	char	*ctrl_file;		//!< Name of control file
	char    *params; 		//!< Solver's comamnd line parameters
	char    keep_files;   	//!< Keep or remove solver files?
	int     use_last_sol;   //!< Use last known solution? (should be in water module)

	//@{
	//! solver parameters
	int             max_it;          //!< Max. number of iterations
	double          r_tol;           //!< Relative tolerance
	double          a_tol;           //!< Absolute tolerance
	ISOL_params		*isol_params;	 //!< optional ISOL parameters
	//@}
} Solver;

// public functions
void solver_init( struct Solver *solver);
void solve_system(struct Solver *solver, LinSys *lin_system);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
