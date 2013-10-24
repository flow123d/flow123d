#ifndef FUNCTORBASE_IMPL_H
#define FUNCTORBASE_IMPL_H

#include "functors.hh"
#include "system/xio.h"

/************************************************ FunctorBase ********************************/
template<class Type>
FunctorBase<Type>::FunctorBase() 
{
  data_.resize(5);   //default resizing of the data vector
}
  
template<class Type>
template<class TType >
FunctorBase<Type>::FunctorBase(FunctorBase<TType>& func) 
{
  unsigned int n = func.n_param();
  data_.resize(n);
  for(unsigned int i=0; i < n; i++) //copying vector
    data_[i] = func.get_param(i);   
}
  
template<class Type>
FunctorBase<Type>::~FunctorBase() {}
  
template<class Type>
void FunctorBase<Type>::set_param(const unsigned int& param_name, const double& param_value)
{
  if(param_name < data_.size())
  {
    data_[param_name] = param_value;
  }
  else
  {
    data_.push_back(param_value);
  }
}
  
template<class Type>
double FunctorBase<Type>::get_param(const unsigned int& param_name)
{
  ASSERT(param_name < data_.size(),"This parameter does not exist.");
  
  return data_[param_name];
}
  
template<class Type>
unsigned int FunctorBase<Type>::n_param()
{ 
  return data_.size();
}
  
  
/************************************************ Functor ************************************/
template<class Type>
Functor<Type>::Functor() : FunctorBase<Type>::FunctorBase()
{}

  
template<class Type>
template<class TType >
Functor<Type>::Functor(Functor<TType>& func) : FunctorBase<Type>::FunctorBase(func)
{}
  
  
template<class Type>
Functor<Type>::~Functor() {}
  
  
/************************************************ FunctorImplicit *****************************/
template<class Type>
FunctorImplicit<Type>::FunctorImplicit() : FunctorBase<Type>::FunctorBase()
{}
  
template<class Type>
template<class TType >
FunctorImplicit<Type>::FunctorImplicit(FunctorImplicit<TType>& func) : FunctorBase<Type>::FunctorBase(func) 
{}
 
  
template<class Type>
FunctorImplicit<Type>::~FunctorImplicit() {}
 

#endif  //FUNCTORBASE_IMPL_H