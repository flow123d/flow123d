#ifndef FUNCTORBASE_IMPL_H
#define FUNCTORBASE_IMPL_H

#include "functors.hh"
#include "system/xio.h"

/************************************************ FunctorBase ********************************/
template<class Type>
FunctorBase<Type>::FunctorBase() 
{}
  
template<class Type>
FunctorBase<Type>::~FunctorBase() {}
  
template<class Type>
void FunctorBase<Type>::set_param(const unsigned int& param_name, const double& param_value)
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
void FunctorBase<Type>::set_param_from_func(FunctorBase<TType>* func)
{
  //DBGMSG("FunctorBase copy constructor.\n");
  unsigned int n = func->n_param();
  //DBGMSG("param_size = %d\n",n);
  param_.resize(n);
  
  for(unsigned int i=0; i < n; i++) //copying vector
  {
    //DBGMSG("get_param = %d\n",i);
    param_[i] = func->param(i);   
  }
}
 
template<class Type>
double FunctorBase<Type>::param(const unsigned int& param_name)
{
  ASSERT(param_name < param_.size(),"Parameter of the functor was not set.");
  
  return param_[param_name];
}
  
template<class Type>
unsigned int FunctorBase<Type>::n_param()
{ 
  return param_.size();
}
  
  
/************************************************ Functor ************************************/
template<class Type>
Functor<Type>::Functor() : FunctorBase<Type>::FunctorBase()
{}
  
  
template<class Type>
Functor<Type>::~Functor() {}
  
  
/************************************************ FunctorImplicit *****************************/
template<class Type>
FunctorImplicit<Type>::FunctorImplicit() : FunctorBase<Type>::FunctorBase()
{}
 
  
template<class Type>
FunctorImplicit<Type>::~FunctorImplicit() {}
 

#endif  //FUNCTORBASE_IMPL_H