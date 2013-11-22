#include "functors_impl.hh"
#include "interpolant_impl.hh"
#include "adaptivesimpson.hh"

#include "system/xio.h"
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
  ASSERT(bound_a != bound_b,"Bounds overlap.");
  ASSERT(bound_a < bound_b,"a must be lower and b must be upper bound.");
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
  ASSERT(size >0, "Size of interpolation table must be positive number!.");
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
  ASSERT(init_size >0, "Maximal size of interpolation table must be positive number!.");
  ASSERT(max_size >= init_size, "Maximal size of interpolation table is smaller than initial size.");
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
  //DBGMSG("Check and reinterpolate.\n");
  bool reinterpolate = false;   //should we remake interpolation?
  if((double)stats.interval_miss_a/stats.total_calls > percentage)
  {
    //DBGMSG("Check and reinterpolate.\n");
    bound_a_ = stats.min;
    reinterpolate = true;
  }
  if((double)stats.interval_miss_b/stats.total_calls > percentage)
  {
    //DBGMSG("Check and reinterpolate.\n");
    bound_b_ = stats.max;
    reinterpolate = true;
  }
  
  if(reinterpolate)
  {
    xprintf(Msg, "Interpolation: Reinterpolating...\n");
    interpolate();
    reset_stat();
  }
}

  
void InterpolantBase::check_all()
{
  ASSERT(checks[Check::functor], "Functor is not set.");
  ASSERT(checks[Check::bound_a], "Left boundary of the interval is not set.");
  ASSERT(checks[Check::bound_b], "Right boundary of the interval is not set.");
  ASSERT(checks[Check::size], "Step is not set.");
  ASSERT(! ((user_tol == 0) && automatic_size), "Tolerance for automatic interpolation is not set.");
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
  ASSERT(n <= n_derivatives,"Not allowed to obtain higher than n-th derivate.");
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
    default: ASSERT(degree > 1, "Higher order interpolation is not available at the moment.");
  }
  */
  
  //temporary vectors for values in the middle of intervals
  std::vector<double> x_vec_half;
  std::vector<double> f_vec_half;
  std::vector<double> df_vec_half;
      
  double tol = 1e-10;   //error computation zero tolerance
  
  if(automatic_size)
  {
    int k = -1;
    create_nodes();     //creates intial set of nodes and values
    while(1) 
    {
      k++;
      //DBGMSG("k = %d\n",k);
      //(this->*interpolate_func)();    //POSSIBLE WAY TO USE MORE KINDS OF INTERPOLATION
      interpolate_p1(); 
      
      if(norm_type == ErrorNorm::max)
        compute_error(tol, x_vec_half, f_vec_half, df_vec_half);
      else
        compute_error(tol, p, norm_type);
      
      DBGMSG("error: %E\n", error_);
    
      //error comparation
      if(user_tol < error_)
      {
        DBGMSG("Interpolating: %d   size: %d \t error: %E\n",k, size_, error_);
        result = 1;   //tolerance has not been satisfied
      }
      else 
      {
        xprintf(Msg,"Interpolation: Size of the table set automaticaly to: %d after %d cycles.\n", size_,k);
        //DBGMSG("Interpolation error estimate is: %f\n", error_);
        result = 0;   //interpolation OK
        break;
      }
      
      //size comparation
      if(size_ < max_size/2)
      {
        if(norm_type == ErrorNorm::max) // if we compute the maximum norm, we can use the computed values
          swap_middle_values(x_vec_half, f_vec_half, df_vec_half); 
        else    //else resize and compute new nodes
        {
          size_ *= 2;
          create_nodes();
        }
      }
      else
      { 
        xprintf(Warn,"Interpolation: User defined tolerance %E has not been satisfied with size of interpolation table %d.\n",user_tol,size_);
        break;
      }
    }
  }
  else 
  {
    //(this->*interpolate_func)();    //POSSIBLE WAY TO USE MORE KINDS OF INTERPOLATION
    create_nodes();     //creates intial set of nodes and values
    interpolate_p1();
    
    if(norm_type == ErrorNorm::max)
      compute_error(tol, x_vec_half, f_vec_half, df_vec_half);
    else
      compute_error(tol, p, norm_type);
      
    result = 0;   //interpolation OK
  }
  xprintf(Msg,"Interpolation: Interpolation error estimate is: %E.\n", error_);
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
*/

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
    /*
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
    //*/
  }

void Interpolant::create_nodes()
{ 
  //setting the step - length of piecewise interpolation intervals
  step = (bound_b_-bound_a_)/size_;
  n_nodes = size_ + 1;          //setting the number of nodes
  
  x_vec.resize(n_nodes);       //nodes
  f_vec.resize(n_nodes);       //function values in the nodes
  
  double temp_x = bound_a_;
  if(interpolate_derivative)
  {
    df_vec.resize(n_nodes);      //function derivates values in the nodes
    //filling the vector x and f
    DiffValue value;
    for(unsigned int i = 0; i < size_; i++)
    {
      //DBGMSG("size: %d \tfill vectors: %d\n",size_,i);
      temp_x = bound_a_ + step*i;
      value = f_diff(temp_x);
      x_vec[i] = temp_x;
      f_vec[i] = value.first;
      df_vec[i] = value.second;
    }
    //finish the interval
    x_vec[size_] = bound_b_;
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
      x_vec[i] = temp_x;
      f_vec[i] = f_val(temp_x);
    }
    //finish the interval
    x_vec[size_] = bound_b_;
    f_vec[size_] = f_val(bound_b_);
  }
}

void Interpolant::swap_middle_values(std::vector<double>& x, std::vector<double>& f, std::vector<double>& df)
{
  double new_size = 2*size_;
        n_nodes = new_size+1;
        step = (bound_b_-bound_a_)/new_size;    //which should be equal also step/2
        
        //we will use now the middle points computed in "compute_error" to construct vector of nodes
        std::vector<double> swap_x_vec;
        std::vector<double> swap_f_vec;
        std::vector<double> swap_df_vec;
        swap_x_vec.reserve(n_nodes);
        swap_f_vec.reserve(n_nodes);

        for(unsigned int i=0; i!=size_; i++)
        {
          swap_x_vec.push_back(x_vec[i]);
          swap_x_vec.push_back(x[i]);
          swap_f_vec.push_back(f_vec[i]);
          swap_f_vec.push_back(f[i]);
        }
        swap_x_vec.push_back(x_vec.back());
        swap_f_vec.push_back(f_vec.back());
        x_vec = swap_x_vec;
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


void Interpolant::compute_error(double tol, std::vector<double>& x, std::vector<double>& f, std::vector<double>& df)
{
  double tot_err = 0;
  
  x.clear();
  f.clear();
  x.reserve(size_);    //middle nodes
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
      //DBGMSG("size: %d \tfill vectors: %d\n",size_,i);
      temp_x = temp_a + step*i;
      value = f_diff(temp_x);
      x.push_back(temp_x);
      f.push_back(value.first);
      df.push_back(value.second);
      
      int_value = diff_p1(temp_x);    //we can call directly the interpolant evaluation, because we know we are in the interval
      //double t = std::abs(value.second - int_value.second) / (std::abs(value.second) + tol);
      tot_err = std::max( tot_err, 
                          std::abs(value.first - int_value.first) / (std::abs(value.first) + tol)
                            + std::abs(value.second - int_value.second) / (std::abs(value.second) + tol)
                        );
      //DBGMSG("x,f,if: \t%f:  %f: %f\n",temp_x, value.first, int_value.first);
      //DBGMSG("temporary tot_error: %f:  %f: %f\n",temp_x,tot_err,t);
    }
  }
  else
  {
    //filling the vector x and f
    for(unsigned int i = 0; i != size_; i++)
    {
      temp_x = temp_a + step*i;
      x.push_back(temp_x);
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
  
  for(unsigned long i = 0; i < p1_vec.size(); i++ )
  { 
    p_err = std::abs( AdaptiveSimpson::AdaptSimpson(*norm,
                                           x_vec[i], 
                                           x_vec[i+1],
                                           simpson_tolerance) );
                   
    //DBGMSG("error on interval<%f,%f>: %f\n",x_vec[i],x_vec[i+1],p_err);
    tot_err += p_err;
  }
  //DBGMSG("tot_err = %f\n",tot_err);
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
    explicit_interpolant(NULL),
    fix_(InterpolantImplicit::no_fix)
  {
  }
  
InterpolantImplicit::~InterpolantImplicit() 
{
  if(explicit_interpolant) delete explicit_interpolant;
  if(func_diff) delete func_diff;
  if(func_diffn) delete func_diffn;
  if(func_u) delete func_u;
}

void InterpolantImplicit::fix_variable(InterpolantImplicit::fix_var fix, const double& val)
{
  fix_ = fix;
  fix_val = val;
  func_u = new InterpolantImplicit::FuncExplicit<double>(*func,fix_,fix_val);
}

double InterpolantImplicit::f_val(double u)
{
  return func_u->operator()(u);
}


void InterpolantImplicit::interpolate_p1()
{
  if(func_u != NULL) 
    delete func_u;
  
  func_u = new FuncExplicit<double>(*func,fix_,fix_val);
}


