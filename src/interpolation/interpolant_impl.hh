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
}


inline double Interpolant::val(const double& i_x)
{
  //increase calls
  stats.total_calls++;
  stats.min = std::min(stats.min, i_x);
  stats.max = std::max(stats.max, i_x);
  
  //DBGMSG("total: %d \t modulo: %d\n",stats.total_calls,stats.total_calls % STATISTIC_CHECK);
  //every STATISTIC_CHECK calls check the statistics
  if((stats.total_calls % STATISTIC_CHECK == 0) && use_statistics) 
    this->check_and_reinterpolate();
  
  //return value
  if(i_x < bound_a_)     //left miss
  {
    stats.interval_miss_a++;
    return f_val(i_x);
  }
  if(i_x > bound_b_)     //right miss
  {
    stats.interval_miss_b++;
    return f_val(i_x);
  }
  else                  //hit the interval
  {
    return val_p1(i_x);
  }
}

inline der Interpolant::diff(const double& i_x)
{
  //increase calls
  stats.total_calls++;
  stats.min = std::min(stats.min, i_x);
  stats.max = std::max(stats.max, i_x);
  
  //every STATISTIC_CHECK calls check the statistics
  if((stats.total_calls % STATISTIC_CHECK == 0) && use_statistics) 
    this->check_and_reinterpolate();
  
  if(i_x < bound_a_)     //left miss
  {
    stats.interval_miss_a++;
    return f_diff(i_x);  //right miss
  }
  if(i_x > bound_b_)
  {
    stats.interval_miss_b++;
    return f_diff(i_x);
  }
  else                  //hit the interval
  {
    return diff_p1(i_x);
  }  
}
  
inline double Interpolant::f_val(const double& i_x)
{
  return func->operator()(i_x);
}
  
inline der Interpolant::f_diff ( const double& i_x )
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
  return floor((i_x - bound_a_) / step);
}

/* //CONSTANT INTERPOLATION
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
*/
  
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



class Interpolant::NormL2 : public Functor<double>
{
public:
  NormL2(Interpolant* interpolant)
  : interpolant(interpolant){}

  virtual double operator()(const double& x)
  {
    return std::pow(interpolant->f_val(x) - interpolant->val(x),2);
  }         
 
private:
  Interpolant* interpolant;
};


class Interpolant::NormW21 : public Functor<double>
{
public:
  NormW21(Interpolant* interpolant)
  : interpolant(interpolant){}
 
  virtual double operator()(const double& x)
  {
    double val = std::pow(interpolant->f_val(x) - interpolant->val(x),2);
    double diff = std::pow(interpolant->f_diff(x).dfdx - interpolant->diff(x).dfdx,2);
    return val+diff;
  }         
  
private:
  Interpolant* interpolant;
};


/********************************** InterpolantImplicit ********************************/

template<template<class> class Func, class Type >
void InterpolantImplicit::set_functor(FunctorImplicit<Type>* func) 
{
  this->func = func;
  func_diff = new Func<B<Type> >(*func);
  func_diffn = new Func<T<Type> >(*func);
  checks[Interpolant::check_functor] = true;
}
  
template<template<class> class Func>
void InterpolantImplicit::set_functor(FunctorImplicit<double>* func) 
{
  this->func = func;
  func_diff = new Func<B<double> >(*func);
  func_diffn = new Func<T<double> >(*func);
  checks[Interpolant::check_functor] = true;
}


  ///class FuncExplicit.
  /** This functor transforms implicit functor with two variables into
   * an explicit functor with only one varible and the other one fixed.
   */
  template<class Type>
  class InterpolantImplicit::FuncExplicit : public Functor<Type>
  {
  public:
    ///Constructor.
    FuncExplicit(){}
  
    //probably not using
    template<class TType>
    FuncExplicit(Functor<TType>& func) : Functor<Type>(func){};
    
    //constructor from templated implicit functor
    template<class TType>
    FuncExplicit(FunctorImplicit<TType>& func_impl, fix_var fix, const double& fix_val)
      : Functor<TType>(func_impl),func_impl(&func_impl), fix_(fix), fix_val(fix_val) {}
  
    virtual Type operator()(const Type& u)
    {
      Type ret;
      if(fix_ == InterpolantImplicit::fix_x)
        ret = func_impl->operator()(fix_val,u);
      if(fix_ == InterpolantImplicit::fix_y)
        ret =  func_impl->operator()(u,fix_val);
      
      return ret;
    }
    
  private:
    FunctorImplicit<Type>* func_impl;
    fix_var fix_;
    double fix_val;
  };