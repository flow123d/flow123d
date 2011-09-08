/*  Sparse matrix subroutines.
    Copyright (C) 2005 Jukka Toivanen
*/

#include "sparse.h"


void CSparseMat::allocate(int rows, int cols, int nnz) {
  irows = rows;
  innz = nnz;
  icols = cols;
  used = 0;

  v = new double[nnz];
  ja = new int[nnz];
  ia = new int[rows+1];

  for (int i=0; i<rows+1; i++)
    ia[i]=0;

  for (int i=0; i<nnz; i++) {
    ja[i]=-1;
    v[i]=0.0;
  }

}  

void CSparseMat::deallocate() {
  delete[] v;
  delete[] ja;
  delete[] ia;
}


CSparseMat::CSparseMat(SuperMatrix *S) {
  if (S->Stype == SLU_NR) {
    NRformat *Store = (NRformat *)S->Store;
    irows = S->nrow;
    icols = S->ncol;
    innz = Store->nnz;
    used = innz;
    v = new double[innz];
    ja = new int[innz];
    ia = new int[irows+1];
      
    for (int i=0; i<innz; i++) {
      v[i] = ((double *)Store->nzval)[i];
      ia[i] = Store->rowptr[i];
      ja[i] = Store->colind[i];
    }
  } else {
    cerr << "Unsupported conversion from SuperMatrix to SparseMat\n";
  }
}



void CSparseMat::Aij(int i, int j, double val, bool add) {
  int k = 0;

  k = ia[i];
  
  do {
    if (k==ia[i+1] || ja[k]>j) {
      if (used>=innz) { cout << "sparse matrix full" << endl; exit(1); }
      for (int h=innz-1; h>k; h--) {
	ja[h]=ja[h-1];
	v[h]=v[h-1];
      }
      for (int h=i+1; h<irows+1; h++)
	ia[h]++;
      ja[k]=j;
      v[k]=val;
      used++;
      break;
    }  

    if (ja[k]==j) {
      if (add) {
	v[k]+=val;
      } else {
	v[k]=val;
      }
      break;
    }
    k=k+1;
  } while(1);
 
}




CSparseMat CSparseMat::operator * (const CSparseMat& B) {
  int rows = irows;
  int cols = B.icols;
  int rownz[cols];
  int nnz = 0;


  for (int i=0; i<rows; i++) {
    for (int j=0; j<cols; j++)
      rownz[j]=0;

    for (int k=ia[i]; k<ia[i+1]; k++)
      for (int j=B.ia[ja[k]]; j<B.ia[ja[k]+1]; j++)
	rownz[B.ja[j]]=1;
    
    for (int j=0; j<cols; j++)
      if (rownz[j]>0) nnz++;
  }
      

  CSparseMat C(rows, cols, nnz); 

  for (int i=0; i<rows; i++) {

    for (int k=ia[i]; k<ia[i+1]; k++)
      for (int j=B.ia[ja[k]]; j<B.ia[ja[k]+1]; j++)
	C.addAij(i, B.ja[j], v[k]*B.v[j]);

  }
  

  return C;

}

/*
  matrix x vector multiplication
*/
double *CSparseMat::operator *(double *vec) {
  double *x = new double[irows];
  for (int i=0; i<irows; i++) x[i] = 0;
  
  for (int i=0; i<irows; i++) {
    for (int j=ia[i]; j<ia[i+1]; j++) {
      x[ja[j]] += v[j]*vec[ja[j]];
    }
  }
  
  return x;
}



void CSparseMat::print() {
  for (int i=0; i<irows; i++) {
    for (int j=ia[i]; j<ia[i+1]; j++) {
      cout << "(" << i <<"," << ja[j] << ") = " << v[j] << endl;
    }
  }
}


void CSparseMat::writeMm(string filename) {
  ofstream f(filename.c_str());
  f << "%%MatrixMarket matrix coordinate real general" << endl;
  f << irows << " " << icols << " " << innz << endl;

  for (int i=0; i<irows; i++) {
    for (int j=ia[i]; j<ia[i+1]; j++) {
      f << i+1 <<" " << ja[j]+1 << " " << v[j] << endl;
    }
  }

  f.close();


}
