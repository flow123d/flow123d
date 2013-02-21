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
 * @brief Auxiliary math functions.
 *
 * Small matrix and vectors should be implemented by armadillo library.
 *
 */

#include <math.h>
#include "global_defs.h"
#include "system/system.hh"
#include "system/math_fce.h"


static double Det2( SmallMtx2 a);
static double Inverse2(SmallMtx2 a,SmallMtx2 b);
static double Inverse3(SmallMtx3 a,SmallMtx3 b);
static double Inverse4(SmallMtx4 a,SmallMtx4 b);

//=============================================================================
// LENGTH OF VECTOR
//=============================================================================
double vector_length( double v[] )
{
	return  sqrt( v[0]*v[0]+v[1]*v[1]+v[2]*v[2] );
}
//=============================================================================
// SCALAR PRODUCT OF VECTORS
//=============================================================================
double scalar_product( double u[], double v[] )
{
	return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
}
//=============================================================================
// NORMALIZE GIVEN VECTOR
//=============================================================================
void normalize_vector( double u[] )
{
	double l;

	if ((l = vector_length( u )) < NUM_ZERO ) {
        xprintf(Warn,"Normalization of nearly zero vector.\n");
	}
	u[0] /= l; u[1] /= l; u[2] /= l;
}
//=============================================================================
// MULTIPLE VECTOR BY REAL NUMBER
//=============================================================================
void scale_vector( double u[], double k )
{
	u[ 0 ] *= k; u[ 1 ] *= k; u[ 2 ] *= k;
}

//=============================================================================
// VECTOR PRODUCT OF VECTORS
//=============================================================================
void vector_product( double u[], double v[], double x[] )
{
	// u, v - inputs
	// x    - output, result

	x[ 0 ] = u[ 1 ] * v[ 2 ] - u[ 2 ] * v[ 1 ];
	x[ 1 ] = u[ 2 ] * v[ 0 ] - u[ 0 ] * v[ 2 ];
	x[ 2 ] = u[ 0 ] * v[ 1 ] - u[ 1 ] * v[ 0 ];
}
//=============================================================================
// VECTOR DIFFERENCE OF VECTORS
//=============================================================================
void vector_difference( double u[], double v[], double x[] )
{
  x[ 0 ] = u[ 0 ] - v[ 0 ];
  x[ 1 ] = u[ 1 ] - v[ 1 ];
  x[ 2 ] = u[ 2 ] - v[ 2 ];
}

/*************************************************
 *  SMALL MATRIX OPERATIONS
 *************************************************/

//=============================================================================
// DETERMINANT OF MATRIX 2x2
//=============================================================================
double Det2( SmallMtx2 a) {
	return SUBDET2(0,1,0,1);
}

/**************************************************
 * determinant 3x3
 * *************************/
double Det3( SmallMtx3 a ) {
	return  a[0][0]*SUBDET2(1,2,1,2)
           -a[0][1]*SUBDET2(1,2,0,2)
		   +a[0][2]*SUBDET2(1,2,0,1);
}

/*************************************************
 * Inverse of general small matrix with given size
 * a - input; b - output; returns determinant
 * ************************************************/
double MatrixInverse(double *a, double *b, int size) {
	switch (size) {
		case 1:
			*b=1/(*a); return (*a);
		case 2:
			return ( Inverse2((SmallMtx2)a, (SmallMtx2)b) );
		case 3:
			return ( Inverse3((SmallMtx3)a, (SmallMtx3)b) );
		case 4:
			return ( Inverse4((SmallMtx4)a, (SmallMtx4)b) );
	}
	return 0.0;
}


/******************************************************************
 * Inverse of 2x2 matrix with determinant
 * a - input; b - output; returns determinant;
 * *****************************************************************/
double Inverse2(SmallMtx2 a,SmallMtx2 b) {
  double Det;

  Det=SUBDET2(0,1,0,1);
  if ( fabs(Det) < NUM_ZERO ) return Det;
  b[0][0]=a[1][1]/Det;
  b[1][0]=-a[1][0]/Det;
  b[0][1]=-a[0][1]/Det;
  b[1][1]=a[0][0]/Det;
  return Det;
}

/******************************************************************
 * Inverse of 3x3 matrix with determinant
 * a - input; b - output; returns determinant;
 * *****************************************************************/

double Inverse3(SmallMtx3 a,SmallMtx3 b) {
	double Det;

	b[0][0]=SUBDET2(1,2,1,2);
	b[1][0]=-SUBDET2(1,2,0,2);
	b[2][0]=SUBDET2(1,2,0,1);

	Det=a[0][0]*b[0][0]+a[0][1]*b[1][0]+a[0][2]*b[2][0];
	if ( fabs(Det)< NUM_ZERO) return Det;
	b[0][0]/=Det; b[1][0]/=Det; b[2][0]/=Det;

	b[0][1]=-SUBDET2(0,2,1,2)/Det;
	b[1][1]=SUBDET2(0,2,0,2)/Det;
	b[2][1]=-SUBDET2(0,2,0,1)/Det;

	b[0][2]=SUBDET2(0,1,1,2)/Det;
	b[1][2]=-SUBDET2(0,1,0,2)/Det;
	b[2][2]=SUBDET2(0,1,0,1)/Det;
	return Det;
}

/******************************************************************
 * Inverse of 4x4 matrix with determinant
 * a - input; b - output; det - determinant; returns determinant
 * *****************************************************************/

double Inverse4(SmallMtx4 a,SmallMtx4 b) {
   double u[6], l[6], Det;
  // 2x2 subdets of rows 1,2
  u[0]=SUBDET2(0,1,0,1);
  u[1]=SUBDET2(0,1,0,2);
  u[2]=SUBDET2(0,1,0,3);
  u[3]=SUBDET2(0,1,1,2);
  u[4]=SUBDET2(0,1,1,3);
  u[5]=SUBDET2(0,1,2,3);

  //1x1 subdets of rows 2,3
  l[0]=SUBDET2(2,3,0,1);
  l[1]=SUBDET2(2,3,0,2);
  l[2]=SUBDET2(2,3,0,3);
  l[3]=SUBDET2(2,3,1,2);
  l[4]=SUBDET2(2,3,1,3);
  l[5]=SUBDET2(2,3,2,3);

  // 3x3 subdets
  b[0][0]=+a[1][1]*l[5]-a[1][2]*l[4]+a[1][3]*l[3];
  b[1][0]=-a[1][0]*l[5]+a[1][2]*l[2]-a[1][3]*l[1];
  b[2][0]=+a[1][0]*l[4]-a[1][1]*l[2]+a[1][3]*l[0];
  b[3][0]=-a[1][0]*l[3]+a[1][1]*l[1]-a[1][2]*l[0];
  // 4x4 det
  Det=a[0][0]*b[0][0]+a[0][1]*b[1][0]+a[0][2]*b[2][0]+a[0][3]*b[3][0];
  if ( fabs(Det)< NUM_ZERO) return Det;
  b[0][0]/=Det; b[1][0]/=Det; b[2][0]/=Det; b[3][0]/=Det;

  b[0][1]=(-a[0][1]*l[5]+a[0][2]*l[4]-a[0][3]*l[3])/Det;
  b[1][1]=(+a[0][0]*l[5]-a[0][2]*l[2]+a[0][3]*l[1])/Det;
  b[2][1]=(-a[0][0]*l[4]+a[0][1]*l[2]-a[0][3]*l[0])/Det;
  b[3][1]=(+a[0][0]*l[3]-a[0][1]*l[1]+a[0][2]*l[0])/Det;

  b[0][2]=(+a[3][1]*u[5]-a[3][2]*u[4]+a[3][3]*u[3])/Det;
  b[1][2]=(-a[3][0]*u[5]+a[3][2]*u[2]-a[3][3]*u[1])/Det;
  b[2][2]=(+a[3][0]*u[4]-a[3][1]*u[2]+a[3][3]*u[0])/Det;
  b[3][2]=(-a[3][0]*u[3]+a[3][1]*u[1]-a[3][2]*u[0])/Det;

  b[0][3]=(-a[2][1]*u[5]+a[2][2]*u[4]-a[2][3]*u[3])/Det;
  b[1][3]=(+a[2][0]*u[5]-a[2][2]*u[2]+a[2][3]*u[1])/Det;
  b[2][3]=(-a[2][0]*u[4]+a[2][1]*u[2]-a[2][3]*u[0])/Det;
  b[3][3]=(+a[2][0]*u[3]-a[2][1]*u[1]+a[2][2]*u[0])/Det;
  return Det;
}

/******************************************************************
 * Print a full matrix stored in a vector size x size
 * *****************************************************************/
void PrintSmallMatrix(double *mtx, int size ) {
	int i,j;
	printf("matrix %d x %d :\n",size,size);
	for (i=0;i<size;i++) {
		for(j=0;j<size;j++) printf("%f ",mtx[i*size+j]);
		printf("\n");
	}
}

/*******************************************************************************
 *  at place inverse of an NxN matrix by Gauss-Jordan elimination without pivotation
 *  this should by faster then via adjung. matrix at least for N>3
 * *****************************************************************************/
void MatInverse(int n,double **a) {
	int row,subrow,col;
	double Div,Mult;

	// forward run - lower triangle and diagonal of the inverse matrix
	for(row=0;row<n-1;row++) {
		// divide the row by diagonal
		Div=a[row][row];a[row][row]=1;
		for(col=0;col<n;col++) a[row][col]/=Div;
		// eliminate the lower part of column "row"
		for(subrow=row+1;subrow<n;subrow++) {
			Mult=a[subrow][row];a[subrow][row]=0;
			for(col=0;col<n;col++) a[subrow][col]-=a[row][col]*Mult;
		}
	}
	// divide the last row
	Div=a[row][row];a[row][row]=1;
	for(col=0;col<n;col++) a[row][col]/=Div;
	//backward run - upper trinagle
	for(;row>0;row--) {
		// eliminate the upper part of column "row"
		for(subrow=row-1;subrow>=0;subrow--) {
			Mult=a[subrow][row];a[subrow][row]=0;
			for(col=0;col<n;col++) a[subrow][col]-=a[row][col]*Mult;
		}
	}
}

//=============================================================================
// solve equation system - gauss; returns 1 - one result 0 - system doesn't have result; 2 - oo results
// !!! THIS FUNCTION IS NOT USED and probably is WRONG !!!
//=============================================================================

int gauss(double *A, double *B,int s,double *R)
{
  int res=1,i,j,k,size;
  double koef,tmp,max;
  double *M;
  size=s+1;
  M = (double *) xmalloc( (size-1) * size * sizeof(double));
  for(i=0;i<s;i++)
  {
    for(j=0;j<s;j++)
      M[i*(size)+j]=A[i*s+j];
    M[i*(size)+j]=B[i];
  }
  for(i=0;i<size-1;i++)
  {
    if(M[i*size+i]==0)
    {
      for(j=i+1;j<size-1;j++)
      {
        if(M[j*size+i]!=0)
        {
          for(k=0;k<size;k++)
          {
             tmp=M[j*size+k];
             M[j*size+k]=M[i*size+k];
             M[i*size+k]=tmp;
          }
          j=size-1;
        }
      }
    }
    max=0;
    for(k=0;k<size;k++)
      if (fabs(max)<fabs(M[i*size+k])) max=fabs(M[i*size+k]);
    if (max!=0)
      for(k=0;k<size;k++)
        M[i*size+k]/=max;
    for(j=i+1;j<size-1;j++)
    {
      if (M[j*size+i]!=0)
      {
        koef=M[i*size+i]/M[j*size+i];
        koef=fabs(koef);
        if((M[i*size+i]>=0 && M[j*size+i]>=0) ||
           (M[i*size+i]<0 && M[j*size+i]<0)) koef*=-1;
        for(k=i;k<size;k++)
          M[j*size+k]=M[j*size+k]*koef + M[i*size+k];
//        M[j*size+i] = 0;
      }
    }
  }
  for(i=size-1-1;i>=0;i--)
  {
    koef=0;
    for(j=i+1;j<size-1;j++)
    {
      koef=koef+R[j]*M[i*size+j];
    }
    if (!(DBL_EQ(M[i*size+i],0)))
//    if (M[i*size+i] != 0)
      R[i]=(M[i*size + size-1]-koef)/M[i*size+i];
    else
    {
      if (koef==0 && M[i*size + size-1]==0)
        return 2;
      else
        return 0;
    }
  }
  xfree(M);
  return res;
}
//=============================================================================
// ANGLE BETWEEN TWO VECTORS
// !!! THIS FUNCTION IS NOT USED !!!
//=============================================================================

double get_vectors_angle( double u[ 3 ], double v[ 3 ] )
{
  double a,b,fi,p;
  a = scalar_product (u, v);
  b = vector_length( u ) * vector_length( v );
  p = a / b;
  if (p > 1) p = 1;
  if (p < -1) p = -1;
  fi = acos( p );
  return fi * 180 / M_PI;
}
//=============================================================================
// MATRIX TIMES MATRIX
// !!! THIS FUNCTION IS NOT USED !!!
//=============================================================================

void matrix_x_matrix(double *A,int ra, int ca,
                     double *B, int rb, int cb, double *X)
{
  int i,j,k;
  ASSERT(!(ca != rb),"Matrix A has different number of columns than matrix B rows in the function matrix_x_matrix()\n");
  for (i = 0; i < ra; i++)
    for (j = 0; j < cb; j++)
    {
      X[ i* cb + j] = 0;
      for (k = 0; k < rb; k++)
        X[ i * cb + j] += A[ i * ca + k] * B[ k * cb + j];
    }
  return;
}
//-----------------------------------------------------------------------------
// vim: set cindent:
