/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    functors_impl.hh
 * @brief   
 */

#ifndef FUNCTORBASE_IMPL_H
#define FUNCTORBASE_IMPL_H

#include "tools/functors.hh"
#include "system/xio.h"

/**************************************** FunctorCommon *****************************************/
template<class Type>
FunctorCommon<Type>::FunctorCommon() {}
  
template<class Type>
FunctorCommon<Type>::~FunctorCommon() {}
  
template<class Type>
void FunctorCommon<Type>::set_param(unsigned int param_name, double param_value)
{
  if(param_name < param_.size())
  {
    param_[param_name] = param_value;
  }
  else
  {
    param_.push_back(param_value);
  }
}

template<class Type> 
template<class TType> 
void FunctorCommon<Type>::set_param_from_func(FunctorCommon<TType>* func)
{
  unsigned int n = func->param_.size();
  //DBGMSG("param_size = %d\n",n);
  param_.resize(n);
  
  for(unsigned int i=0; i < n; i++) //copying vector
  {
    //DBGMSG("get_param = %d\n",i);
    param_[i] = func->param(i);   
  }
}
 
template<class Type>
double FunctorCommon<Type>::param(unsigned int param_name)
{
  OLD_ASSERT(param_name < param_.size(),"Parameter of the functor was not set.");
  
  return param_[param_name];
}

  
/************************************************ FunctorBase ************************************/
template<class Type>
FunctorBase<Type>::FunctorBase() : FunctorCommon<Type>::FunctorCommon() {}
  
template<class Type>
FunctorBase<Type>::~FunctorBase() {}
  
  
/************************************************ IFunctorBase ***********************************/
template<class Type>
IFunctorBase<Type>::IFunctorBase() : FunctorCommon<Type>::FunctorCommon() {}
  
template<class Type>
IFunctorBase<Type>::~IFunctorBase() {}
 

#endif  //FUNCTORBASE_IMPL_H
