#include <cmath>
#include "mesh/ngh/include/mathfce.h"


bool mathfce::IsZero(double x)
{
  if (fabs(x) < epsilon)
    return true;
  return false;
}

bool mathfce::IsEqual(double x, double y)
{
  if (fabs(x - y) < epsilon)
    return true;
  return false;
}


double Determinant3( double a[ 3 ][ 3 ] )
{
  double rc;
  rc =   a[ 0 ][ 0 ] * a[ 1 ][ 1 ] * a[ 2 ][ 2 ]
       + a[ 1 ][ 0 ] * a[ 2 ][ 1 ] * a[ 0 ][ 2 ]
       + a[ 2 ][ 0 ] * a[ 0 ][ 1 ] * a[ 1 ][ 2 ]
       - a[ 0 ][ 2 ] * a[ 1 ][ 1 ] * a[ 2 ][ 0 ]
       - a[ 1 ][ 2 ] * a[ 2 ][ 1 ] * a[ 0 ][ 0 ]
       - a[ 2 ][ 2 ] * a[ 0 ][ 1 ] * a[ 1 ][ 0 ];
  return rc;
}

