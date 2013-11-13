#include "interpolant.hh"
#include "system/xio.h"

#include <cmath>

/********************************** InterpolantBase ********************************/

inline double InterpolantBase::error()
{
  return error_;
}

inline InterpolantBase::EvalStatistics InterpolantBase::statistics() const
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
  
  checks[Check::functor] = true;
}


template<template<class> class Func, class Type >
void Interpolant::set_functor(Func<Type>* func) 
{
  this->func = func;
  func_diff = new Func<B<Type> >();
  func_diffn = new Func<T<Type> >();
  
  func_diff->set_param_from_func(func);
  func_diffn->set_param_from_func(func);
  
  checks[Check::functor] = true;
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
    
    switch(extrapolation)
    {
      case Extrapolation::constant: 
        return f_vec[0];
      case Extrapolation::linear:
        return f_vec[0] + p1_vec[0]*(x-x_vec[0]);
      case Extrapolation::functor:
      default:
        return f_val(x);      //otherwise compute original functor
    }
  }
  else if(x > bound_b_)     //right miss
  {
    stats.interval_miss_b++;
    stats.max = std::max(stats.max, x);
    
    switch(extrapolation)
    {
      case Extrapolation::constant: 
        return f_vec[size_];
      case Extrapolation::linear:
        return f_vec[size_-1] + p1_vec[size_-1]*(x-x_vec[size_-1]);
      case Extrapolation::functor:
      default:
        return f_val(x);      //otherwise compute original functor
    }
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
    switch(extrapolation)
    {
      case Extrapolation::constant: 
        return DiffValue(f_vec[0],
                         df_vec[0]);
      case Extrapolation::linear:
        return DiffValue(f_vec[0] + p1_vec[0]*(x-x_vec[0]), 
                         df_vec[0] + p1d_vec[0]*(x-x_vec[0]));
      case Extrapolation::functor:
      default:
        return f_diff(x);      //otherwise compute original functor
    }
  }
  else if(x > bound_b_)      //right miss
  {
    stats.interval_miss_b++;
    stats.max = std::max(stats.max, x);
    switch(extrapolation)
    {
      case Extrapolation::constant: 
        return DiffValue(f_vec[size_], 
                         df_vec[size_]);
      case Extrapolation::linear:
        return DiffValue(f_vec[size_-1] + p1_vec[size_-1]*(x-x_vec[size_-1]),
                         df_vec[size_-1] + p1d_vec[size_-1]*(x-x_vec[size_-1]) );
      case Extrapolation::functor:
      default:
        return f_diff(x);      //otherwise compute original functor
    }
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



/** Functor class that computes the argument of the integral \f$ (f(x)-i(x))^2 \f$ in the norm \f$ \|f-i\|_{L_2}\f$. 
   * It is used as input functor to integration.
   */
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

  /** Functor class that computes the argument of the integral \f$ (f(x)-i(x))^2 + (f'(x)-i'(x))^2 \f$ in the norm \f$ \|f-i\|_{W^1_2} \f$. 
   * It is used as input functor to integration.
   */
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
        
  /** Functor class that computes the argument of the integral \f$ (f(x)-i(x))^p + (f'(x)-i'(x))^p \f$ in the norm \f$ \|f-i\|_{W^1_p} \f$. 
   * It is used as input functor to integration.
   */
class Interpolant::NormWp1 : public FunctorBase<double>
{
public:
  NormWp1(Interpolant* interpolant, double p)
  : interpolant(interpolant), p_(p) {}
 
  inline double p() {return p_;}
 
  virtual double operator()(double x)
  {
    double val = std::pow(interpolant->f_val(x) - interpolant->val(x),p_);
    double diff = std::pow(interpolant->f_diff(x).second - interpolant->diff(x).second,p_);
    return val+diff;
  }         
  
private:
  Interpolant* interpolant;
  double p_;
};

  /** Functor class that computes the argument of the integral \f$ (f(x)-i(x))^p + (f'(x)-i'(x))^p \f$ in the norm \f$ \|f-i\|_{W^1_p} \f$. 
   * It is used as input functor to integration.
   */
class Interpolant::NormFunc : public FunctorBase<double>
{
public:
  NormFunc(Interpolant* interpolant, double p)
  : interpolant(interpolant), p_(p) {}
 
  inline double p() {return p_;}
 
  virtual double operator()(double x)
  {
    double val = std::pow(interpolant->f_val(x),p_);
    double diff = std::pow(interpolant->f_diff(x).second,p_);
    return val+diff;
  }         
  
private:
  Interpolant* interpolant;
  double p_;
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
  
  checks[Check::functor] = true;
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