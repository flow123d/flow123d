/* Solution of systems of linear equations
   Jukka Toivanen, Jan Stebel
*/

#ifndef _SOLUTION_H__
#define _SOLUTION_H__

#include <fstream>
#include <string>
#include <iostream>
#include <vector>

#include "slu_ddefs.h"
#include "umfpack.h"
#include "ad.h"


using namespace std;


/**
  Solves the equation A x = b using library SuperLU.
  @param mat   Sparse matrix.
  @param x     Array of right hand size (on return - solution vector).
  @param trans Flag indicating solution with transposed matrix.
*/
int superSolve   (CSparseMat &mat, double *x, bool trans);

/**
  Solves the equation A x = b usinglibrary UMFPACK.
  @param mat   Sparse matrix.
  @param x     Array of right hand size (on return - solution vector).
  @param trans Flag indicating solution with transposed matrix.
*/
int UMFPACKSolve(CSparseMat &mat, double *x, bool trans);

/* unused
  Solves dr/dq delta = r and sets q = q - relax * delta.

int solveState(CVec &r, CVec &q, double relax, int low, int ncomp);
int solveAdjoint(CVec &r, CVar &J, CVec &q, int low, int ncomp, double *p);
void saveMatrix(CVec &r, CVec &q, int low, int ncomp, char *file);
*/
#endif
