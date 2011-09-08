/* Automatic differentiation subroutines - Forward mode.
   Jukka Toivanen, Jan Stebel
*/

#ifndef _AD_H__
#define _AD_H__

#include <iostream>
#include "sparse.h"

/** Default size of array of derivatives. */
#define MAX_AD 40

typedef struct
{
  int ind;
  double fd, gd;
} cder;

/**
  \brief Basic class for automatic differentiation.
  
  The class manipulates double-typed variables together with their first
  derivatives with respect to independent variables.
*/
class CVar
{
 private:
  static int count;  /**< Actual number of CVar variables. */
  int        ind,    /**< Index of the variable. */
             nder,   /**< Number of stored derivatives. */
	     asize;  /**< Maximal number of derivatives. */
  double*    der;    /**< Array of derivatives. */
  int*       dermap; /**< Array of variables with respect to which
                          the derivatives are taken. */
  double     val;    /**< The value of the variable. */
  
  int        fi, gi;
  
  inline cder getFirstCder(const CVar&);
  inline cder getNextCder(const CVar&);

 public:
  /** Default constructor. */
  CVar();
  
  /**
    Constructor initializing the variable with a given value.
    @param value Initial value.
  */
  CVar(double value);
  
  /**
    Constructor initializing the variable with a given value
    and a given size of the array of derivatives.
    @param value Initial value.
    @param size  Size of the array of derivatives.
  */
  CVar(double value, int size);
  
  /** Copy constructor. */
  CVar(const CVar& var);
  
  /** Destructor. */
  ~CVar();
  
  /**
    Allocates the arrays and index.
    @param size Array size.
  */
  void   init(int size);
  
  /**
    Sets derivative with respect to index-th variable to value.
    @param index Derivative index.
    @param value Derivative value.
  */
  void   setDer(int index, double value);
  
  /** Returns Derivative with respect to index-th variable.
      @param index Derivative index. */
  double getDer(int index) const;
  
  /**
    Returns vector of derivatives in a specified range.
    @param low   Starting index of derivatives.
    @param ncomp Number of derivatives to be returned.
                 Ending index of derivatives is low+ncomp.
    @param der   Array of derivatives. Must be allocated before calling
                 the routine.
  */
  void   getDer(int low, int ncomp, double* der) const;
  
  /**
    Returns the derivative with respect to a given variable.
    @param var Differentiating variable.
  */
  double getDerOfVar(const CVar& var) const;
  
  /**
    Returns the derivative with respect to variable given by index.
    @param Index of differentiating variable.
  */
  double getDerOfVar(int ind) const;
  
  /**
    Returns the index of variable in the stored array of derivatives.
    @param ind Position within the array of indices of differentiating
               variables.
  */
  int getDerVar(int ind) const;
  
  /** Returns index of the actual variable. */
  int getInd() const;
  
  /** Returns number of derivatives of the actual variable. */
  int getNder() const;
  
  /** Returns the value of the actual variable. */
  inline double getVal() const { return val; };
  
  /**
    Sets the actual variable independent so that derivatives w.r.t this
    variable are calculated
  */
  void setIndependent();
  
  /** Clear the array of derivatives. */
  void resetDer();
  
  CVar  operator + (double);
  CVar  operator + (const CVar&);
  CVar& operator +=(const CVar&);
  CVar  operator - (double);
  CVar  operator - (const CVar&);
  CVar  operator * (double);
  CVar  operator * (const CVar&);
  CVar  operator / (const CVar&);
  CVar  operator - ();
  CVar& operator = (double);
  CVar& operator = (const CVar&);
  
  friend double sgn(double x);
  friend CVar abs(CVar var);
  friend CVar xabsx(CVar var);
  friend CVar sqrt(CVar var);
  friend CVar sin(CVar var);
  friend CVar cos(CVar var);
  friend CVar tanh(CVar var);
  friend CVar pow(CVar base, double exponent);
  friend CVar pow(double base, CVar exponent);
};



class CVec {
 private:
  CVar* vec;
  int asize;
 public:
  CVec();
  CVec(int size);
  CVec(const CVec &);
  ~CVec();
  double norm2();
  double min();
  double max();
  void resize(int);
  void setIndependent();
  void resetDer();
  int size();
  CVar& operator() (int loc);
  CVar operator | (const CVec&);
  CVec& operator = (const CVec&);
  CVec& operator = (const double&);
  void getVal(double *val) const;
  void getJac(int low, int ncomp, CSparseMat &jac) const;
  void getSubJac(int M, int m, int low, int ncomp, CSparseMat &jac) const;
};


template <typename T> class CMat
{
 private:
  T *mat;
  int nr, nc;
  void del();
 public:
  CMat();
  CMat(int rows, int cols);
  ~CMat();
  void resize(int rows, int cols);
  void resetDer();
  int getNrows() { return nr; };
  int getNcols() { return nc; };
  T norm2();
  T& operator () (int row, int col);
  CMat<T>& operator = (const CMat<T>&);
};

template <typename T> CMat<T>::CMat()
{
  nr=0;
  nc=0;
  mat = new T[nr*nc];
}

template <typename T> CMat<T>::CMat(int rows, int cols)
{
  nr = rows;
  nc = cols;
  mat = new T[nr*nc];
}

template <typename T> CMat<T>::~CMat()
{
  del();
}

template <typename T> void CMat<T>::del()
{
  delete[] mat;
}

template <typename T> void CMat<T>::resize(int rows, int cols)
{
  if (nr!=rows || nc!=cols) {
    del();
    nr = rows;
    nc = cols;
    mat = new T[nr*nc];
  }
}

template <typename T> T CMat<T>::norm2()
{
  T norm = 0;
  int i,j;
  for (i=0;i<nr;i++)
    for (j=0;j<nc;j++)
      norm += mat[i*nc+j]*mat[i*nc+j];
      
  return sqrt(norm);
}

template <typename T> T& CMat<T>::operator() (int row, int col)
{
  return mat[row*nc+col];
}

template <typename T> CMat<T>& CMat<T>::operator = (const CMat<T>& param)
{
  resize(param.nr, param.nc);
  for (int i=0; i<nr*nc; i++)
    mat[i]=param.mat[i];

  return *this;
}



template <typename T> void mat_inverse_3d(CMat<T> &J, CMat<T> &IJ, T &dJ, T &IdJ)
{
  dJ = J(0,0)*J(1,1)*J(2,2) + J(1,0)*J(2,1)*J(0,2) + J(2,0)*J(0,1)*J(1,2)
      -J(0,0)*J(2,1)*J(1,2) - J(1,0)*J(0,1)*J(2,2) - J(2,0)*J(1,1)*J(0,2);
  IdJ=1.0;
  IdJ=IdJ/dJ;
  
  IJ(0,0)=IdJ*(J(1,1)*J(2,2)-J(2,1)*J(1,2));
  IJ(0,1)=IdJ*(J(0,2)*J(2,1)-J(2,2)*J(0,1));
  IJ(0,2)=IdJ*(J(0,1)*J(1,2)-J(1,1)*J(0,2));
  IJ(1,0)=IdJ*(J(1,2)*J(2,0)-J(2,2)*J(1,0));
  IJ(1,1)=IdJ*(J(0,0)*J(2,2)-J(2,0)*J(0,2));
  IJ(1,2)=IdJ*(J(0,2)*J(1,0)-J(1,2)*J(0,0));
  IJ(2,0)=IdJ*(J(1,0)*J(2,1)-J(2,0)*J(1,1));
  IJ(2,1)=IdJ*(J(0,1)*J(2,0)-J(2,1)*J(0,0));
  IJ(2,2)=IdJ*(J(0,0)*J(1,1)-J(1,0)*J(0,1));
}

// nasledujici funkce spocita pseudo-inverzi matice J typu 2x3, dJ = sqrt(det(J^T*J))
template <typename T> void mat_pseudoinverse_2x3(CMat<T> &J, CMat<T> &IJ, T &dJ, T &IdJ)
{
  dJ  = (J(0,0)*J(0,0)+J(0,1)*J(0,1)+J(0,2)*J(0,2))*(J(1,0)*J(1,0)+J(1,1)*J(1,1)+J(1,2)*J(1,2))
       -(J(0,0)*J(1,0)+J(0,1)*J(1,1)+J(0,2)*J(1,2))*(J(0,0)*J(1,0)+J(0,1)*J(1,1)+J(0,2)*J(1,2));
  IdJ = (T)1.0/dJ;
  
  IJ(0,0) = IdJ*(J(0,0)*J(1,1)*J(1,1) + J(0,0)*J(1,2)*J(1,2) - J(1,0)*J(0,1)*J(1,1) - J(1,0)*J(0,2)*J(1,2));
  IJ(0,1) = IdJ*(J(1,0)*J(0,1)*J(0,1) + J(1,0)*J(0,2)*J(0,2) - J(0,0)*J(0,1)*J(1,1) - J(0,0)*J(0,2)*J(1,2));
  IJ(1,0) = IdJ*(J(0,1)*J(1,0)*J(1,0) + J(0,1)*J(1,2)*J(1,2) - J(0,0)*J(1,0)*J(1,1) - J(0,2)*J(1,1)*J(1,2));
  IJ(1,1) = IdJ*(J(0,0)*J(0,0)*J(1,1) + J(1,1)*J(0,2)*J(0,2) - J(0,0)*J(1,0)*J(0,1) - J(0,1)*J(0,2)*J(1,2));
  IJ(2,0) = IdJ*(J(1,0)*J(1,0)*J(0,2) + J(1,1)*J(1,1)*J(0,2) - J(0,0)*J(1,0)*J(1,2) - J(0,1)*J(1,1)*J(1,2));
  IJ(2,1) = IdJ*(J(0,0)*J(0,0)*J(1,2) + J(0,1)*J(0,1)*J(1,2) - J(0,0)*J(1,0)*J(0,2) - J(0,1)*J(1,1)*J(0,2));
  
  dJ = sqrt(dJ);
}

template <typename T> void mat_inverse_2d(CMat<T> &J, CMat<T> &IJ, T &dJ, T &IdJ)
{
  dJ  = J(0,0)*J(1,1)-J(1,0)*J(0,1);
  IdJ = 1.0;
  IdJ = IdJ/dJ;
  IJ(0,0)=IdJ*J(1,1);
  IJ(0,1)=-IdJ*J(0,1);
  IJ(1,0)=-IdJ*J(1,0);
  IJ(1,1)=IdJ*J(0,0);
}

template <typename T> void mat_pseudoinverse_1x3(CMat<T> &J, CMat<T> &IJ, T &dJ, T &IdJ)
{
  dJ  = J(0,0)*J(0,0)+J(0,1)*J(0,1)+J(0,2)*J(0,2);
  IdJ = 1.0;
  IdJ = IdJ/dJ;
  IJ(0,0)=IdJ*J(0,0);
  IJ(1,0)=IdJ*J(0,1);
  IJ(2,0)=IdJ*J(0,2);
  dJ = sqrt(dJ);
}





#endif
