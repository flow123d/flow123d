/*  Sparse matrix subroutines.
    Copyright (C) 2005 Jukka Toivanen
*/

#ifndef _SPARSE_H__
#define _SPARSE_H__


#include <iostream>
#include <fstream>
#include "slu_ddefs.h"

using namespace std;

/**
  Class for handling sparse matrices.
  A[i, ja[j]] = v[ia[i]+j]
*/
class CSparseMat
{
 private:
  int    *ja,   /**< Array of column indices for the given entry in v. */
         *ia;   /**< Array of pointers to first entries in given row.  */
  double *v;    /**< Array of non-zero entries.                        */
  int    used,  /**< Number of used nonzero values. */
         irows, /**< Number of rows. */
         icols, /**< Number of columns. */
         innz;  /**< Number of max. nonzero values. */
  
  /**
    Set matrix coefficient.
    @param i   Row index.
    @param j   Column index.
    @param val Value of coefficient
    @param add If true, val be added to current coefficient value.
  */
  void Aij(int i, int j, double val, bool add);
  
  /**
    Allocate matrix arrays.
    @param rows Number of rows.
    @param cols Number of columns.
    @param nnz  Number of nonzero elements.
  */
  void allocate(int rows, int cols, int nnz);
  
  /** Destroy allocated arrays. */
  void deallocate();
 public:
  /** Inplicit constructor. */
  CSparseMat() { allocate(0,0,0);  }
  
  /**
    Standard constructor.
    @param rows Number of rows.
    @param cols Number of columns.
    @param nnz  Number of nonzero elements.
  */
  CSparseMat(int rows, int cols, int nnz) { allocate(rows, cols, nnz); }
  
  /** Constructor for converting from SuperMatrix format (used by SuperLU). */
  CSparseMat(SuperMatrix*);
  
  /** Destructor. */
  ~CSparseMat() { deallocate(); }
  
  /**
    Add value to matrix coefficient.
    @param i   Row index.
    @param j   Column index.
    @param val Value to be added.
  */
  void addAij(int i, int j, double val) { Aij(i, j, val, true); };
  
  /**
    Set matrix coefficient.
    @param i   Row index.
    @param j   Column index.
    @param val Value to be set.
  */
  void setAij(int i, int j, double val) { Aij(i, j, val, false); };
  
  /** Getter for ia. */
  int* getIa() { return ia; }
  
  /* Getter for ja. */
  int* getJa() { return ja; }
  
  /** Return number of columns. */
  int getNcols() { return icols; }
  
  /** Return max. number of nonzero values. */
  int getNnz() { return innz; }
  
  /** Return number of rows. */
  int getNrows() { return irows; }
  
  /**
    Resize matrix.
    @param rows New number of rows.
    @param cols New number of columns.
    @param nnz  New number of nonzero values.
  */
  void resize(int rows, int cols, int nnz) { deallocate(); allocate(rows, cols, nnz); }
  
  /** Print matrix data. */
  void print();
  
  /** Getter for v. */
  double* getV() { return v; }
  
  /** Write in matrix-market format. */
  void writeMm(string filename);
  
  /** Matrix multiplication. */
  CSparseMat operator * (const CSparseMat&);
  
  /** Multiplication by scalar. */
  double *operator * (double *);
  
  friend class CVec;
};

#endif
