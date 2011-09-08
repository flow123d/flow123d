/* Solution of systems of linear equations
   Jukka Toivanen, Jan Stebel
*/

#include <math.h>
#include "solution.h"

int nprocs = 1;

int superSolve(CSparseMat &mat, double *x, bool trans) {
     SuperMatrix A, L, U, B, X;
     int      *perm_r; /* row permutations from partial pivoting */
     int      *perm_c; /* column permutation vector */
     int      info, m, n, nnz;
     superlu_options_t options;
     SuperLUStat_t stat;
     int            *etree;
     char           equed[1];
     void           *work = NULL;
     int lwork;
     double rpg, rcond;
     double         *ferr, *berr, *rhsx, *xact;
     double         *R, *C;
     mem_usage_t    mem_usage;


     lwork = 0;

     m = n = mat.getNrows();

     nnz = mat.getNnz();

     /* Create matrix A in the format expected by SuperLU. */
     dCreate_CompRow_Matrix(&A, m, n, nnz, mat.getV(), mat.getJa(), mat.getIa(), SLU_NR, SLU_D, SLU_GE);

     dCreate_Dense_Matrix(&B, m, 1, x, m, SLU_DN, SLU_D, SLU_GE);

     if ( !(rhsx = doubleMalloc(m)) ) ABORT("Malloc fails for rhsx[].");
     dCreate_Dense_Matrix(&X, m, 1, rhsx, m, SLU_DN, SLU_D, SLU_GE);


    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
     if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) )
         ABORT("SUPERLU_MALLOC fails for R[].");
     if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
         ABORT("SUPERLU_MALLOC fails for C[].");
     if ( !(ferr = (double *) SUPERLU_MALLOC(sizeof(double))) )
         ABORT("SUPERLU_MALLOC fails for ferr[].");
     if ( !(berr = (double *) SUPERLU_MALLOC(sizeof(double))) )
         ABORT("SUPERLU_MALLOC fails for berr[].");
     xact = doubleMalloc(n);

     /* Set the default input options. */
     set_default_options(&options);

     if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
     if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");

     

     if (trans)
       options.Trans = TRANS;

     /* Initialize the statistics variables. */
     StatInit(&stat);


     dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
            &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
            &mem_usage, &stat, &info);
     
     if (info!=0) { cerr << "Linear system is not solvable, info " << info << endl; return 1; }


     for (int i=0; i<n; i++)
       x[i]=rhsx[i];


     /* De-allocate storage */

     SUPERLU_FREE (rhsx);
     SUPERLU_FREE (xact);
     SUPERLU_FREE (etree);
     SUPERLU_FREE (perm_r);
     SUPERLU_FREE (perm_c);
     SUPERLU_FREE (R);
     SUPERLU_FREE (C);
     SUPERLU_FREE (ferr);
     SUPERLU_FREE (berr);
     Destroy_SuperMatrix_Store(&A);
     Destroy_SuperMatrix_Store(&B);
     Destroy_SuperMatrix_Store(&X);
     if ( lwork >= 0 ) {
         Destroy_SuperNode_Matrix(&L);
         Destroy_CompCol_Matrix(&U);
     }

     StatFree(&stat);
     
     return 0;
}


int UMFPACKSolve(CSparseMat &mat, double *x, bool trans)
{
  double *null = (double *) NULL, *Ax, *sol;
  int n, nz, *Ap, *Ai;
  void *Symbolic, *Numeric ;
  
  n  = mat.getNrows();
  nz = mat.getNnz();
  Ap = mat.getIa();
  Ai = mat.getJa();
  Ax = mat.getV();
  
  sol = new double[n];
  
  (void) umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, null, null) ;
  (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
  umfpack_di_free_symbolic (&Symbolic) ;
  if (trans)
  {
    (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, sol, x, Numeric, null, null) ;
  }
  else
  {
    (void) umfpack_di_solve (UMFPACK_At, Ap, Ai, Ax, sol, x, Numeric, null, null) ;
  }
  umfpack_di_free_numeric (&Numeric) ;
  
  for (int i=0; i<n; i++)
  {
    x[i] = sol[i];
  }
  delete[] sol;
  
  return 0;
}


int solveState(CVec &r, CVec &q, double relax, int low, int ncomp) {
    CSparseMat sp;
//    cout << "compute r jacobian...\n";
    r.getJac(q(low).getInd(), ncomp, sp);
//    cout << "ok\n";
    int n = sp.getNrows();
    
    double *x = new double[n];

    for (int i=0; i<n; i++)
      x[i] = r(i).getVal();

    if (UMFPACKSolve(sp, x, false) != 0) return 1;

    for (int i=0; i<n; i++)
      q(i) = q(i)-x[i]*relax;

    delete[] x;
    return 0;
}


int solveAdjoint(CVec &r, CVar &J, CVec &q, int low, int ncomp, double *p) {
  CSparseMat sp;
  
  J.getDer(q(low).getInd(), ncomp, p);

  r.getJac(q(low).getInd(), ncomp, sp);
  
  int res = superSolve(sp, p, true);
  
  return res;

}


void saveMatrix(CVec &r, CVec &q, int low, int ncomp, char *file)
{
  ofstream f(file);
#ifndef AD_REV
  for (int i=0; i<r.size(); i++)
  {
    for (int j=0; j<r(i).getNder(); j++)
    {
      if (r(i).getDer(j) != 0 && r(i).getDerVar(j) >= q(low).getInd() && r(i).getDerVar(j) < q(low).getInd()+ncomp)
      {
        f << i << " " << r(i).getDerVar(j)-q(low).getInd() << " " << r(i).getDer(j) << endl;
      }
    }
  }
#endif
  f.close();
}





