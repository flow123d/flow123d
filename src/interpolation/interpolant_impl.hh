#ifndef INTERPOLANT_IMPL_H
#define INTERPOLANT_IMPL_H

#include <cmath>

#include "interpolant.hh"
#include "system/xio.h"

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
Interpolant::Interpolant(Func<Type>* func, bool interpolate_derivative) 
  : func(func), 
    interpolate_derivative(interpolate_derivative)
{
  func_diff = new Func<B<Type> >();
  func_diffn = new Func<T<Type> >();
  
  func_diff->set_param_from_func(func);
  func_diffn->set_param_from_func(func);
  
  checks[Check::functor] = true;
}


template<template<class> class Func, class Type >
void Interpolant::set_functor(Func<Type>* func, bool interpolate_derivative) 
{
  this->interpolate_derivative = interpolate_derivative;
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
    //DBGMSG("test\n");
    //std::cout << "test a" << std::endl;
    stats.interval_miss_a++;
    stats.min = std::min(stats.min, x);
    
    switch(extrapolation)
    {
      case Extrapolation::constant: 
        return f_vec[0];
      case Extrapolation::linear:
        return (f_vec[1] - f_vec[0])*(x/step-a_div_step) + f_vec[0]; //f_vec[0] + p1_vec[0]*(x-x_vec[0]);
      case Extrapolation::functor:
      default:
        return f_val(x);      //otherwise compute original functor
    }
  }
  else if(x > bound_b_)     //right miss
  {
    //DBGMSG("test\n");
    //std::cout << "test b" << std::endl;
    stats.interval_miss_b++;
    stats.max = std::max(stats.max, x);
    
    switch(extrapolation)
    {
      case Extrapolation::constant: 
        return f_vec[size_];
      case Extrapolation::linear:
        return (f_vec[size_] - f_vec[size_-1])*(x-bound_b_+step)/step + f_vec[size_-1];//f_vec[size_-1] + p1_vec[size_-1]*(x-x_vec[size_-1]);
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

inline double Interpolant::val_test(double x)
{
  /*
  //increase calls
  stats.total_calls++;
 
  if(x < bound_a_)
  {
    //std::cout << "test a" << std::endl;
    return 0;
  }
  else if(x > bound_b_)     //right miss
  {
    //std::cout << "test b" << std::endl;
    return 0;
  }
  else  //*/
     return val_p1(x);  
}

inline DiffValue Interpolant::diff(double x)
{
  ASSERT(interpolate_derivative, "Derivative is not interpolated. Flag must be switched true in constructor (or set_functor).");
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
        return DiffValue((f_vec[1] - f_vec[0])*(x/step-a_div_step) + f_vec[0], //f_vec[0] + p1_vec[0]*(x-x_vec[0]), 
                          (df_vec[1] - df_vec[0])*(x/step-a_div_step) + df_vec[0]);//df_vec[0] + p1d_vec[0]*(x-x_vec[0]));
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
        return DiffValue((f_vec[size_] - f_vec[size_-1])*(x-bound_b_+step)/step + f_vec[size_-1], //f_vec[size_-1] + p1_vec[size_-1]*(x-x_vec[size_-1]),
                          (df_vec[size_] - df_vec[size_-1])*(x-bound_b_+step)/step + df_vec[size_-1]);//df_vec[size_-1] + p1d_vec[size_-1]*(x-x_vec[size_-1]) );
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

/*
inline unsigned int Interpolant::find_interval(double x)
{
  //counts in which interval x is (the last node before x)
  return floor((x - bound_a_) / step);
}
*/

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
  /*
[       OK ] InterpolantTest.InterpolantEval_val (184 ms)
[       OK ] InterpolantTest.InterpolantEval_val_p1 (192 ms)
[       OK ] InterpolantTest.InterpolantEval_val_test (190 ms)
  
  When we used precomputed coeficients of P1 interpolation.
   */
  //unsigned int i = floor(x / step - a_div_step);
  //return p1_vec[i]*(x-x_vec[i]) + f_vec[i];
  
  /*
[       OK ] InterpolantTest.InterpolantEval_val (162 ms)
[       OK ] InterpolantTest.InterpolantEval_val_p1 (163 ms)
[       OK ] InterpolantTest.InterpolantEval_val_test (162 ms)

These results are satisfying:
 + the function val_p1 itself is by 14% faster than before
 + we use only one array to keep the function values
 + all three functions take the same time which means
   that the collecting statistics and going through several IF commands
   does not affect the evaluation time
  */
  double x_step = x / step;
  unsigned int i = floor(x_step - a_div_step);
  
  return (f_vec[i+1] - f_vec[i])*(x_step-i-a_div_step) + f_vec[i];
}


inline DiffValue Interpolant::diff_p1(double x)
{
  DiffValue result;
  //unsigned int i = find_interval(x);
  //result.first = p1_vec[i]*(x-x_vec[i]) + f_vec[i];
  //result.second = p1d_vec[i]*(x-x_vec[i]) + df_vec[i];
  
  double x_step = x / step;
  unsigned int i = floor(x_step - a_div_step);
  
  result.first =  (f_vec[i+1] - f_vec[i])*(x_step-i-a_div_step) + f_vec[i];
  result.second =  (df_vec[i+1] - df_vec[i])*(x_step-i-a_div_step) + df_vec[i];
  
  return result;
}


/** Functor class that computes the argument of the integral \f$ (f(x)-i(x))^p + (f'(x)-i'(x))^p \f$ in the norm \f$ \|f-i\|_{W^1_p} \f$. 
   * It is used as input functor to integration.
   */
class Interpolant::FuncError_lp : public FunctorBase<double>
{
public:
  FuncError_lp(Interpolant* interpolant, double p, double tol)
  : interpolant(interpolant), p_(p), tol_(tol) {}
 
  inline double p() {return p_;}
  inline double tol() {return tol_;}
 
  virtual double operator()(double x)
  {
    double f_val = interpolant->f_val(x);
    double a = std::abs(f_val - interpolant->val(x)) / (std::abs(f_val) + tol_);
             
    return std::pow(a, p_);
  }         
  
private:
  Interpolant* interpolant;
  double p_;
  double tol_;
};
       
  /** Functor class that computes the argument of the integral \f$ (f(x)-i(x))^p + (f'(x)-i'(x))^p \f$ in the norm \f$ \|f-i\|_{W^1_p} \f$. 
   * It is used as input functor to integration.
   */
class Interpolant::FuncError_wp1 : public FunctorBase<double>
{
public:
  FuncError_wp1(Interpolant* interpolant, double p, double tol)
  : interpolant(interpolant), p_(p), tol_(tol) {}
 
  inline double p() {return p_;}
  inline double tol() {return tol_;}
 
  virtual double operator()(double x)
  {
    DiffValue f = interpolant->f_diff(x);
    DiffValue g = interpolant->diff(x);
    double a = std::abs(f.first - g.first) / (std::abs(f.first) + tol_)
             + std::abs(f.second - g.second) / (std::abs(f.second) + tol_);
             
    return std::pow(a, p_);
  }         
  
private:
  Interpolant* interpolant;
  double p_;
  double tol_;
};



/********************************** InterpolantImplicit ********************************/

template<template<class> class Func, class Type >
void InterpolantImplicit::set_functor(Func<Type>* func, bool interpolate_derivative) 
{
  this->func = func;
  func_diff = new Func<B<Type> >();
  func_diffn = new Func<T<Type> >();
  
  func_diff->set_param_from_func(func);
  func_diffn->set_param_from_func(func);
  
  this->interpolate_derivative = interpolate_derivative;
  
  checks[Check::functor] = true;
}

template<template<class> class Func, class Type >
InterpolantImplicit::InterpolantImplicit(Func<Type>* func, bool interpolate_derivative)
  : func(func), 
    interpolate_derivative(interpolate_derivative),
    explicit_interpolant(NULL),
    fix_(IFixVariable::no_fix)
{
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
    FuncExplicit(IFunctorBase<TType>& func_impl, IFixVariable::Type fix, double fix_val)
      : func_impl(&func_impl), fix_(fix), fix_val(fix_val) {}
    
    Type operator()(Type u)
    {
      Type ret;
      if(fix_ == IFixVariable::fix_x)
        ret = func_impl->operator()(fix_val,u);
      if(fix_ == IFixVariable::fix_y)
        ret =  func_impl->operator()(u,fix_val);
      
      return ret;
    }
    
  private:
    IFunctorBase<Type>* func_impl;
    IFixVariable::Type fix_;
    double fix_val;
  };
  
  
inline double InterpolantImplicit::val(double u)
{
  return explicit_interpolant->val(u);
}

inline DiffValue InterpolantImplicit::diff(double u)
{
  return explicit_interpolant->diff(u);
}

#endif //INTERPOLATION_IMPL_H