#include "interpolant.hh"

/********************************** Interpolant ********************************/  

template<template<class> class Func, class Type >
void Interpolant::set_functor(Functor<Type>* func) 
{
  this->func = func;
  func_diff = new Func<B<Type> >(*func);
  func_diffn = new Func<T<Type> >(*func);
  checks[Interpolant::check_functor] = true;
}
  
template<template<class> class Func>
void Interpolant::set_functor(Functor<double>* func) 
{
  this->func = func;
  func_diff = new Func<B<double> >(*func);
  func_diffn = new Func<T<double> >(*func);
  checks[Interpolant::check_functor] = true;
  
  DBGMSG("Functor copy constructor: param(2)=%f\n",func_diff->get_param(2));
}


inline double Interpolant::val(const double& i_x)
{
  if(i_x < bound_a)
  {
    interval_miss_a++;
    return f_val(i_x);
  }
  if(i_x > bound_b)
  {
    interval_miss_b++;
    return f_val(i_x);
  }
  else
  {
    interval_hits++;
    return (this->*val_)(i_x);
  }
}

inline der Interpolant::diff(const double& i_x)
{
  if(i_x < bound_a)
  {
    interval_miss_a++;
   return f_diff(i_x);
  }
  if(i_x > bound_b)
  {
    interval_miss_b++;
    return f_diff(i_x);
  }
  else
  {
    interval_hits++;
    return (this->*diff_)(i_x);
  }  
}
  
inline double Interpolant::f_val(const double& i_x)
{
  return func->operator()(i_x);
}
  
inline  der Interpolant::f_diff ( const double& i_x )
{
  B<double> x(i_x);   // Initialize arguments
  //Func func;        // Instantiate functor
  B<double> f(func_diff->operator()(x)); // Evaluate function and record DAG
  f.diff(0,1);        // Differentiate
 
  der d;
  d.f = f.x();        // Value of function
  d.dfdx = x.d(0);    // Value of df/dx
  return d;           // Return function value
}

inline unsigned int Interpolant::find_interval(const double& i_x)
{
  //counts in which interval x is (the last node before x)
  return floor((i_x - bound_a) / step);
}

inline double Interpolant::val_p0(const double& i_x)
{
  return f_vec[find_interval(i_x)];
}

inline der Interpolant::diff_p0(const double& i_x)
{
  der result;
  unsigned int i = find_interval(i_x);
  result.f = f_vec[i];
  result.dfdx = df_vec[i];
  return result;
}
  
inline double Interpolant::val_p1(const double& i_x)
{
  unsigned int i = find_interval(i_x);
  return p1_vec[i]*(i_x-x_vec[i]) + f_vec[i];
}
  
inline der Interpolant::diff_p1(const double& i_x)
{
  der result;
  unsigned int i = find_interval(i_x);
  result.f = p1_vec[i]*(i_x-x_vec[i]) + f_vec[i];
  result.dfdx = p1d_vec[i]*(i_x-x_vec[i]) + df_vec[i];
  return result;
}




/********************************** InterpolantImplicit ********************************/    

