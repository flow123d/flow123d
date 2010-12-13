//---------------------------------------------------------------------------
//    $Id: tensor.cc 14038 2006-10-23 02:46:34Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <tensor.h>
#include <cmath>


DEAL_II_NAMESPACE_OPEN


// storage for static variables
template <int dim>
const unsigned int Tensor<1,dim>::dimension;

template <int rank, int dim>
const unsigned int Tensor<rank,dim>::dimension;


template class Tensor<1, 1>;
template class Tensor<1, 2>;
template class Tensor<1, 3>;
template class Tensor<2, 1>;
template class Tensor<2, 2>;
template class Tensor<2, 3>;
template class Tensor<2, 4>;
//template class Tensor<3, 1>;
//template class Tensor<3, 2>;
//template class Tensor<3, 3>;
//template class Tensor<4, 1>;
//template class Tensor<4, 2>;
//template class Tensor<4, 3>;

DEAL_II_NAMESPACE_CLOSE
