#include "interpolant.hh"
#include "system/xio.h"

#include <cmath>

/********************************** InterpolantBase ********************************/

inline double InterpolantBase::error()
{
  return error_;
}

inline InterpolantBase::eval_statistics InterpolantBase::statistics() const
{
  return stats;
}

inline double InterpolantBase::bound_a() const
{
  return bound_a_;
}

inline double InterpolantBase::bound_b() const
{
  return bound_b_;
}

inline unsigned int InterpolantBase::size() const
{
  return size_;
}


/********************************** Interpolant ********************************/

template<template<class> class Func, class Type >
Interpolant::Interpolant(Func<Type>* func) 
{
  this->func = func;
  func_diff = new Func<B<Type> >();
  func_diffn = new Func<T<Type> >();
  
  func_diff->set_param_from_func(func);
  func_diffn->set_param_from_func(func);
  
  checks[Interpolant::check_functor] = true;
}


template<template<class> class Func, class Type >
void Interpolant::set_functor(Func<Type>* func) 
{
  this->func = func;
  func_diff = new Func<B<Type> >();
  func_diffn = new Func<T<Type> >();
  
  func_diff->set_param_from_func(func);
  func_diffn->set_param_from_func(func);
  
  checks[Interpolant::check_functor] = true;
}


inline double Interpolant::val(double x)
{
  //increase calls
  stats.total_calls++;
  
  //uncomment if we want to do the check automatically
  //every STATISTIC_CHECK calls check the statistics
  //if((stats.total_calls % STATISTIC_CHECK == 0) && use_statistics) 
  //  this->check_and_reinterpolate();
  
  //return value
  if(x < bound_a_)     //left miss
  {
    stats.interval_miss_a++;
    stats.min = std::min(stats.min, x);
    return f_val(x);
  }
  if(x > bound_b_)     //right miss
  {
    stats.interval_miss_b++;
    stats.max = std::max(stats.max, x);
    return f_val(x);
  }
  else                  //hit the interval
  {
    return val_p1(x);
  }
}

inline DiffValue Interpolant::diff(double x)
{
  //increase calls
  stats.total_calls++;
  
  //uncomment if we want to do the check automatically
  //every STATISTIC_CHECK calls check the statistics
  //if((stats.total_calls % STATISTIC_CHECK == 0) && use_statistics) 
  //  this->check_and_reinterpolate();
  
  if(x < bound_a_)     //left miss
  {
    stats.interval_miss_a++;
    stats.min = std::min(stats.min, x);
    return f_diff(x);  //right miss
  }
  if(x > bound_b_)
  {
    stats.interval_miss_b++;
    stats.max = std::max(stats.max, x);
    return f_diff(x);
  }
  else                  //hit the interval
  {
    return diff_p1(x);
  }  
}
  
inline double Interpolant::f_val(double x)
{
  return func->operator()(x);
}
  
inline DiffValue Interpolant::f_diff (double x)
{
  B<double> xx(x);   // Initialize arguments
  //Func func;        // Instantiate functor
  B<double> f(func_diff->operator()(xx)); // Evaluate function and record DAG
  f.diff(0,1);        // Differentiate
 
  DiffValue d;
  d.first = f.x();    // Value of function
  d.second = xx.d(0);  // Value of df/dx
  return d;           // Return function value
}

inline unsigned int Interpolant::find_interval(double x)
{
  //counts in which interval x is (the last node before x)
  return floor((x - bound_a_) / step);
}

/* //CONSTANT INTERPOLATION
inline double Interpolant::val_p0(double x)
{
  return f_vec[find_interval(x)];
}

inline DiffValue Interpolant::diff_p0(double x)
{
  DiffValue result;
  unsigned int i = find_interval(x);
  result.f = f_vec[i];
  result.dfdx = df_vec[i];
  return result;
}
*/
  
inline double Interpolant::val_p1(double x)
{
  unsigned int i = find_interval(x);
  return p1_vec[i]*(x-x_vec[i]) + f_vec[i];
}
  
inline DiffValue Interpolant::diff_p1(double x)
{
  DiffValue result;
  unsigned int i = find_interval(x);
  result.first = p1_vec[i]*(x-x_vec[i]) + f_vec[i];
  result.second = p1d_vec[i]*(x-x_vec[i]) + df_vec[i];
  return result;
}



class Interpolant::NormL2 : public FunctorBase<double>
{
public:
  NormL2(Interpolant* interpolant)
  : interpolant(interpolant){}

  virtual double operator()(double x)
  {
    return std::pow(interpolant->f_val(x) - interpolant->val(x),2);
  }         
 
private:
  Interpolant* interpolant;
};


class Interpolant::NormW21 : public FunctorBase<double>
{
public:
  NormW21(Interpolant* interpolant)
  : interpolant(interpolant){}
 
  virtual double operator()(double x)
  {
    double val = std::pow(interpolant->f_val(x) - interpolant->val(x),2);
    double diff = std::pow(interpolant->f_diff(x).second - interpolant->diff(x).second,2);
    return val+diff;
  }         
  
private:
  Interpolant* interpolant;
};


/********************************** InterpolantImplicit ********************************/

template<template<class> class Func, class Type >
void InterpolantImplicit::set_functor(Func<Type>* func) 
{
  this->func = func;
  func_diff = new Func<B<Type> >();
  func_diffn = new Func<T<Type> >();
  
  func_diff->set_param_from_func(func);
  func_diffn->set_param_from_func(func);
  
  checks[Interpolant::check_functor] = true;
}

  ///class FuncExplicit.
  /** This functor transforms implicit functor with two variables into
   * an explicit functor with only one variable and the other one fixed.
   */
  template<class Type>
  class InterpolantImplicit::FuncExplicit : public FunctorBase<Type>
  {
  public:
    ///Constructor.
    FuncExplicit(){}
    
    //constructor from templated implicit functor
    template<class TType>
    FuncExplicit(IFunctorBase<TType>& func_impl, fix_var fix, double fix_val)
      : func_impl(&func_impl), fix_(fix), fix_val(fix_val) {}
    
    virtual Type operator()(Type u)
    {
      Type ret;
      if(fix_ == InterpolantImplicit::fix_x)
        ret = func_impl->operator()(fix_val,u);
      if(fix_ == InterpolantImplicit::fix_y)
        ret =  func_impl->operator()(u,fix_val);
      
      return ret;
    }
    
  private:
    IFunctorBase<Type>* func_impl;
    fix_var fix_;
    double fix_val;
  };