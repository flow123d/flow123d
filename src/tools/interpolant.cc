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
 * @file    interpolant.cc
 * @brief   
 */

#include "adaptivesimpson.hh"
#include "functors_impl.hh"
#include "interpolant_impl.hh"
#include "system/logger.hh"
//#include "system/sys_profiler.hh"

#include <limits>
#include <iomanip>

/********************************** InterpolantBase ******************************/
const unsigned int InterpolantBase::n_derivatives = 10;
const unsigned int InterpolantBase::default_max_size = 10*1000;
const double InterpolantBase::simpson_tolerance = 1e-10;

InterpolantBase::InterpolantBase()
  : bound_a_(0),
    bound_b_(0),
    max_size(default_max_size),
    automatic_size(false),
    norm_type(ErrorNorm::max),
    error_(-1.0),
    extrapolation(Extrapolation::functor)
    //use_statistics(true)
{
  reset_stat();
  checks.resize(4,false);
}

InterpolantBase::~InterpolantBase() 
{}
  

void InterpolantBase::set_interval(double bound_a, double bound_b)
{
  OLD_ASSERT(bound_a != bound_b,"Bounds overlap.");
  OLD_ASSERT(bound_a < bound_b,"a must be lower and b must be upper bound.");
  this->bound_a_ = bound_a;
  checks[Check::bound_a] = true;
  this->bound_b_ = bound_b;
  checks[Check::bound_b] = true;
  
  //used to be use to collect the minimum and maximum "x" through all the evaluations (not only outside interval)
  //are given oposite to be able to response to all calls
  //stats.min = std::numeric_limits<double>::max();
  //stats.max = -std::numeric_limits<double>::max();
  stats.min = bound_a;
  stats.max = bound_b;
}
  
void InterpolantBase::set_size(unsigned int size)
{
  OLD_ASSERT(size >0, "Size of interpolation table must be positive number!.");
  this->size_ = size;
  checks[Check::size] = true;
  
  automatic_size = false;
}

void InterpolantBase::set_norm(ErrorNorm::Type norm_type, double p)
{
  this->norm_type = norm_type;
  this->p = p;
}


void InterpolantBase::set_size_automatic(double user_tol, unsigned int init_size, unsigned int max_size)
{
  OLD_ASSERT(init_size >0, "Maximal size of interpolation table must be positive number!.");
  OLD_ASSERT(max_size >= init_size, "Maximal size of interpolation table is smaller than initial size.");
  this->user_tol = user_tol;
  this->max_size = max_size;
  size_ = init_size;
  automatic_size = true;
  
  checks[Check::size] = true;
}

void InterpolantBase::set_extrapolation(Extrapolation::Type extrapolation)
{
  this->extrapolation = extrapolation;
}

void InterpolantBase::reset_stat()
{
  stats.interval_miss_a = 0;
  stats.interval_miss_b = 0;
  stats.total_calls = 0;
  //stats.min = std::numeric_limits<double>::max();
  //stats.max = -std::numeric_limits<double>::max();
  stats.min = bound_a_;
  stats.max = bound_b_;
  //switch statistics on
  //use_statistics = true;
}

void InterpolantBase::check_stats_and_reinterpolate(double percentage)
{
  bool reinterpolate = false;   //should we remake interpolation?
  if((double)stats.interval_miss_a/stats.total_calls > percentage)
  {
    bound_a_ = stats.min;
    reinterpolate = true;
  }
  if((double)stats.interval_miss_b/stats.total_calls > percentage)
  {
    bound_b_ = stats.max;
    reinterpolate = true;
  }
  
  if(reinterpolate)
  {
    MessageOut() << "Interpolation: Reinterpolating...\n";
    interpolate();
    reset_stat();
  }
}

  
void InterpolantBase::check_all()
{
  OLD_ASSERT(checks[Check::functor], "Functor is not set.");
  OLD_ASSERT(checks[Check::bound_a], "Left boundary of the interval is not set.");
  OLD_ASSERT(checks[Check::bound_b], "Right boundary of the interval is not set.");
  OLD_ASSERT(checks[Check::size], "Step is not set.");
  OLD_ASSERT(! ((user_tol == 0) && automatic_size), "Tolerance for automatic interpolation is not set.");
}

long InterpolantBase::fact(long n)
{
  if (n > 1)
    return n*fact(n-1);
  else
    return 1;
}
  
/********************************** Interpolant ********************************/  
  
Interpolant::Interpolant()
  : InterpolantBase(),
    func(NULL),
    func_diff(NULL),
    func_diffn(NULL),
    interpolate_derivative(false)
  {
  }
  
Interpolant::~Interpolant() 
{
  if(func_diff) delete func_diff;
  if(func_diffn) delete func_diffn;
}

  
double Interpolant::f_diffn(double x,unsigned int n)
{
  OLD_ASSERT(n <= n_derivatives,"Not allowed to obtain higher than n-th derivate.");
  T<double> xx,f;

  xx = x;
  xx[1] = 1;
  f = func_diffn->operator()(xx);
  f.eval(n_derivatives);
      
  for(unsigned int i = 0; i<n_derivatives; i++)    //goes through the field of Taylor's coeficients
  {                           //and divide them by the factorial to get derivates
    f[i] = f[i] * this->fact(i);
  }
    
  return f[n];                //returns n-th derivate
}

int Interpolant::interpolate()
{   
  check_all();  //checks if all needed parameters were set
  unsigned int result = 10;

  //switch off statistics due to error computation
  //is switched on automatically by reset_stat()
  //use_statistics = false;

  //POSSIBLE WAY TO USE MORE KINDS OF INTERPOLATION
  //selecting the interpolation 
  /*
  void (Interpolant::*interpolate_func)(void);
  switch(degree)
  {
    case 0: interpolate_func = &Interpolant::interpolate_p0; break; 
    case 1: interpolate_func = &Interpolant::interpolate_p1; break; 
    default: OLD_ASSERT(degree > 1, "Higher order interpolation is not available at the moment.");
  }
  */
  
  //temporary vectors for values in the middle of intervals
  std::vector<double> f_vec_half;
  std::vector<double> df_vec_half;
      
  double tol = 1e-10;   //error computation zero tolerance
  
  if(automatic_size)
  {
    int k = -1;
    while(1) 
    {
      k++;
      compute_values();     //computes function (and possibly derivate) values at nodes
      
      if(norm_type == ErrorNorm::max)
        compute_error(tol, f_vec_half, df_vec_half);
      else
        compute_error(tol, p, norm_type);
      
      DebugOut().fmt("error: {}\n", error_);
    
      //error comparation
      if(user_tol < error_)
      {
        DebugOut().fmt("Interpolating: {}   size: {} \t error: {}\n",k, size_, error_);
        result = 1;   //tolerance has not been satisfied
      }
      else 
      {
        MessageOut().fmt("Interpolation: Size of the table set automaticaly to: {} after {} cycles.\n", size_, k);
        //DebugOut().fmt("Interpolation error estimate is: {}\n", error_);
        result = 0;   //interpolation OK
        break;
      }
      
      //size comparation
      if(size_ < max_size/2)
      {
        if(norm_type == ErrorNorm::max) // if we compute the maximum norm, we can use the computed values
          swap_middle_values(f_vec_half, df_vec_half); 
        else    //else resize and compute new nodes
        {
          size_ *= 2;
          compute_values();     //computes function (and possibly derivate) values at nodes
        }
      }
      else
      { 
        WarningOut().fmt("Interpolation: User defined tolerance {} has not been satisfied with size of interpolation table {}).\n",
        		user_tol, size_);
        break;
      }
    }
  }
  else 
  {
    compute_values();     //computes function (and possibly derivate) values at nodes

    if(norm_type == ErrorNorm::max)
      compute_error(tol, f_vec_half, df_vec_half);
    else
      compute_error(tol, p, norm_type);
      
    result = 0;   //interpolation OK
  }
  MessageOut().fmt("Interpolation: Interpolation error estimate is: {}.\n", error_);
  reset_stat();
  return result;
}

/*
  // CONSTANT INTERPOLATION
  void Interpolant::interpolate_p0()
  {
    create_nodes();

    Functor<double>* norm = new NormW21(this);
    compute_error(norm);      //sets error_
    delete norm;
  }

  // LINEAR INTERPOLATION WITH SAVED COEFICIENTS
  void Interpolant::interpolate_p1()
  {
    p1_vec.resize(size_);    //linear coeficients

    if(interpolate_derivative)
    {
      double delta = 0;
      p1d_vec.resize(size_);    //linear coeficients for derivative
      for(unsigned int i = 0; i != p1_vec.size(); i++)
      {
        delta = x_vec[i+1] - x_vec[i];
        p1_vec[i] = (f_vec[i+1] - f_vec[i]) / delta;
        p1d_vec[i] = (df_vec[i+1] - df_vec[i]) / delta;
      }
    }
    else
    {
      for(unsigned int i = 0; i != p1_vec.size(); i++)
      {
        p1_vec[i] = (f_vec[i+1] - f_vec[i]) / (x_vec[i+1] - x_vec[i]);
      }
    }

    //Writes the interpolation table.

    unsigned int p = 10,
                 pp = 5;
    for(unsigned int i=0; i != size_-1; i++)
    {
      std::cout << "x: " << setw(p) << x_vec[i] << setw(pp) << "f:" << setw(p) << f_vec[i] << setw(pp) << "p1:" << setw(p) << p1_vec[i];
      if(interpolate_derivative)
        std::cout << setw(pp) << "df:" << setw(p) << df_vec[i] << setw(pp) << "p1d:" << setw(p) << p1d_vec[i] << std::endl;
      else
        std::cout << std::endl;
    }
      std::cout << "x: " << setw(p) << x_vec[x_vec.size()-1] << setw(pp) << "f:" << setw(p) << f_vec[f_vec.size()-1] << setw(pp) << "p1:" << setw(p) << "-";
      if(interpolate_derivative)
        std::cout << setw(pp) << "df:" << setw(p) << df_vec[df_vec.size()-1] << setw(pp) <<  "p1d:" << setw(p) << "-" << std::endl;
      else
        std::cout << std::endl;
  }
//*/

void Interpolant::compute_values()
{ 
  //setting the step - length of piecewise interpolation intervals
  step = (bound_b_-bound_a_)/size_;
  a_div_step = bound_a_ / step;
  n_nodes = size_ + 1;          //setting the number of nodes
  
  f_vec.resize(n_nodes);       //function values in the nodes
  
  double temp_x = bound_a_;
  if(interpolate_derivative)
  {
    df_vec.resize(n_nodes);      //function derivates values in the nodes
    //filling the vector x and f
    DiffValue value;
    for(unsigned int i = 0; i < size_; i++)
    {
      temp_x = bound_a_ + step*i;
      value = f_diff(temp_x);
      f_vec[i] = value.first;
      df_vec[i] = value.second;
    }
    //finish the interval
    value = f_diff(bound_b_);
    f_vec[size_] = value.first;
    df_vec[size_] = value.second;
  }
  else
  {
    //filling the vector x and f
    for(unsigned int i = 0; i < size_; i++)
    {
      temp_x = bound_a_ + step*i;
      f_vec[i] = f_val(temp_x);
    }
    //finish the interval
    f_vec[size_] = f_val(bound_b_);
  }
}

void Interpolant::swap_middle_values(std::vector<double>& f, std::vector<double>& df)
{
  double new_size = 2*size_;
        n_nodes = new_size+1;
        step = (bound_b_-bound_a_)/new_size;    //which should be equal also step/2
        a_div_step = bound_a_ / step;
        
        //we will use now the middle points computed in "compute_error" to construct vector of nodes
        std::vector<double> swap_f_vec;
        std::vector<double> swap_df_vec;
        swap_f_vec.reserve(n_nodes);

        for(unsigned int i=0; i!=size_; i++)
        {
          swap_f_vec.push_back(f_vec[i]);
          swap_f_vec.push_back(f[i]);
        }
        swap_f_vec.push_back(f_vec.back());
        f_vec = swap_f_vec;
        
        if(interpolate_derivative)
        {
          swap_df_vec.reserve(n_nodes);
          for(unsigned int i=0; i!=size_; i++)
          {
            swap_df_vec.push_back(df_vec[i]);
            swap_df_vec.push_back(df[i]);
          }
          swap_df_vec.push_back(df_vec.back());
          df_vec = swap_df_vec;
        } 
        size_ = new_size; 
}


void Interpolant::compute_error(double tol, std::vector<double>& f, std::vector<double>& df)
{
  double tot_err = 0;
  
  f.clear();
  f.reserve(size_);    //function values in the middle nodes

  double temp_a = bound_a_ + step/2;
  double temp_x = temp_a;
  if(interpolate_derivative)
  {
    df.clear();
    df.reserve(size_); //function derivates values in the middle nodes
    //filling the vector x and f
    DiffValue value;
    DiffValue int_value;
    for(unsigned int i = 0; i != size_; i++)
    {
      temp_x = temp_a + step*i;
      value = f_diff(temp_x);
      f.push_back(value.first);
      df.push_back(value.second);
      
      int_value = diff_p1(temp_x);    //we can call directly the interpolant evaluation, because we know we are in the interval
      //double t = std::abs(value.second - int_value.second) / (std::abs(value.second) + tol);
      tot_err = std::max( tot_err, 
                          std::abs(value.first - int_value.first) / (std::abs(value.first) + tol)
                            + std::abs(value.second - int_value.second) / (std::abs(value.second) + tol)
                        );
    }
  }
  else
  {
    //filling the vector x and f
    for(unsigned int i = 0; i != size_; i++)
    {
      temp_x = temp_a + step*i;
      f.push_back(f_val(temp_x));
      tot_err = std::max( tot_err, 
                          std::abs(f.back() - val_p1(temp_x)) / (std::abs(f.back()) + tol)
                        );
    }
  }
  
  error_ = tot_err;
  
  /* PRIORITY QUEUE FOR ADAPTIVE INTERVAL DIVIDING
  for(unsigned long i = 0; i < g->get_count(); i++ )
  {
    LP_Norm norm(f,g->get_polynomial(i),2);
    p_err.i = i;

    //absolute polynomial error
    p_err.err = sqrt(Interpolation::AdaptiveSimpson::AdaptSimpson( norm,
                                              g->get_polynomial(i)->get_a(),
                                              g->get_polynomial(i)->get_b(),
                                              SIMPSON_TOLERANCE) );

    //increase the absolute total error
    tot_err += p_err.err;

    //writes absolute error on a single polynomial
    if(DEB) std::cout << "\t abs. p_err=" << p_err.err;

    //p_err convertion absolute -> relative (p_err/(xi+1 - xi))
    p_err.err /= (g->get_polynomial(i)->get_b()-g->get_polynomial(i)->get_a());

    //writes relative error on a single polynomial
    if(DEB) std::cout << "\t rel. p_err=" << p_err.err << std::endl;

    pq.push(p_err); //puts in priority queue
  }

  // writes relative and absolute total error
  //std::cout << "\trelative err = "
  //  << tot_err/(g->GetA()-g->GetB())
  //  << "\tabsolute err = " << tot_err << std::endl;
  //*/
}

void Interpolant::compute_error(double tol, double p, ErrorNorm::Type norm_type)
{
  double exponent = 2;  //default exponent
  FunctorBase<double>* norm;
  
  switch(norm_type)
  {     
    case ErrorNorm::w21:    
        norm = new FuncError_wp1(this, 2, tol);
        break;
    case ErrorNorm::wp1:    
        norm = new FuncError_wp1(this, p, tol);
        exponent = p;
        break;
    case ErrorNorm::lp:    
        norm = new FuncError_lp(this, p, tol);
        exponent = p;
        break;
    case ErrorNorm::l2:    
    default:
        norm = new FuncError_lp(this, 2, tol);                    
  }
  

  double  tot_err = 0, // total absolute error on <a,b>
          p_err = 0;   // absolute error on x[i]-x[i+1]
  double x_node;
  
  for(unsigned long i = 0; i < size_; i++ )
  { 
    x_node = bound_a_ + i*step;
    p_err = std::abs( AdaptiveSimpson::AdaptSimpson(*norm,
                                           x_node, 
                                           x_node+step,
                                           simpson_tolerance) );
                   
    tot_err += p_err;
  }
  tot_err /= (bound_b_ - bound_a_);       //relative to the interval length  
  tot_err = std::pow(tot_err, 1.0/exponent);     //p-th root
  error_ = tot_err;
}
  
  
  
/********************************** InterpolantImplicit ********************************/  
  
InterpolantImplicit::InterpolantImplicit()
  : InterpolantBase(),
    func_u(NULL),
    func(NULL),
    func_diff(NULL),
    func_diffn(NULL),
    interpolate_derivative(false),
    explicit_interpolant(NULL),
    fix_(IFixVariable::no_fix)
  {
  }
  
InterpolantImplicit::~InterpolantImplicit() 
{
  if(explicit_interpolant) delete explicit_interpolant;
  if(func_diff) delete func_diff;
  if(func_diffn) delete func_diffn;
  if(func_u) delete func_u;
}

void InterpolantImplicit::fix_variable(IFixVariable::Type fix, double val)
{
  fix_ = fix;
  fix_val = val;
  func_u = new InterpolantImplicit::FuncExplicit<double>(*func,fix_,fix_val);
}

double InterpolantImplicit::f_val(double u)
{
  return func_u->operator()(u);
}


int InterpolantImplicit::interpolate()
{
  //BUG: FuncExplicit cannot be copied in interpolant constructor with its members !!!
  /*
  OLD_ASSERT(fix_ != IFixVariable::no_fix,"Cannot do interpolation. No varible was fixed.");
  check_all();
  DebugOut().fmt("seg {}\n",func_u);
  explicit_interpolant = new Interpolant(func_u, interpolate_derivative);
  //explicit_interpolant->set_functor(func_u, interpolate_derivative);
  explicit_interpolant->set_interval(bound_a_,bound_b_);
  if(automatic_size)
  {
    explicit_interpolant->set_size_automatic(user_tol, size_, max_size);
  }
  else
  {
    explicit_interpolant->set_size(size_);
  }

  explicit_interpolant->interpolate();
  error_ = explicit_interpolant->error();
  */
  return 0;
}

void InterpolantImplicit::interpolate_p1()
{
  if(func_u != NULL) 
    delete func_u;
  
  func_u = new FuncExplicit<double>(*func,fix_,fix_val);
}


