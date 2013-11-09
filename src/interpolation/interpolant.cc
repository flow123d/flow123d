#include "functors_impl.hh"
#include "interpolant_impl.hh"
#include "adaptivesimpson.hh"

#include "system/xio.h"
//#include "system/sys_profiler.hh"

#include <limits>

/********************************** InterpolantBase ******************************/
const unsigned int InterpolantBase::n_derivatives = 10;

InterpolantBase::InterpolantBase()
  : bound_a_(0),
    bound_b_(0),
    max_size(MAX_SIZE),
    automatic_step(false),
    error_(-1.0),
    extrapolation(InterpolantBase::functor),
    use_statistics(true)
{
  reset_stat();
  checks.resize(4);
}
  
InterpolantBase::~InterpolantBase() 
{}
  

void InterpolantBase::set_interval(double bound_a, double bound_b)
{
  ASSERT(bound_a != bound_b,"Bounds overlap.");
  ASSERT(bound_a < bound_b,"a must be lower and b must be upper bound.");
  this->bound_a_ = bound_a;
  checks[check_a] = true;
  this->bound_b_ = bound_b;
  checks[check_b] = true;
  
  //are given oposite be able to response to all calls
  //stats.min = std::numeric_limits<double>::max();
  //stats.max = -std::numeric_limits<double>::max();
  stats.min = bound_a;
  stats.max = bound_b;
}
  
void InterpolantBase::set_size(unsigned int size)
{
  ASSERT(size >0, "Size of interpolation table must be positive number!.");
  this->size_ = size;
  checks[Interpolant::check_size] = true;
  
  automatic_step = false;
}

void InterpolantBase::set_size_automatic(double user_tol, unsigned int init_size, unsigned int max_size)
{
  ASSERT(init_size >0, "Maximal size of interpolation table must be positive number!.");
  ASSERT(max_size >= init_size, "Maximal size of interpolation table is smaller than initial size.");
  this->user_tol = user_tol;
  this->max_size = max_size;
  size_ = init_size;
  automatic_step = true;
}

void InterpolantBase::set_extrapolation(InterpolantBase::ExtrapolationType extrapolation)
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
  //switch statistics on
  stats.min = bound_a_;
  stats.max = bound_b_;
  use_statistics = true;
}

void InterpolantBase::check_and_reinterpolate()
{
  //DBGMSG("Check and reinterpolate.\n");
  bool reinterpolate = false;   //should we remake interpolation?
  if((double)stats.interval_miss_a/stats.total_calls > MISS_PERCENTAGE)
  {
    DBGMSG("Check and reinterpolate.\n");
    bound_a_ = stats.min;
    reinterpolate = true;
  }
  if((double)stats.interval_miss_b/stats.total_calls > MISS_PERCENTAGE)
  {
    //DBGMSG("Check and reinterpolate.\n");
    bound_b_ = stats.max;
    reinterpolate = true;
  }
  
  if(reinterpolate)
  {
    DBGMSG("Reinterpolating...\n");
    interpolate();
    reset_stat();
  }
}

  
void InterpolantBase::check_all()
{
  ASSERT(checks[Interpolant::check_functor], "Functor is not set.");
  ASSERT(checks[Interpolant::check_a], "Left boundary of the interval is not set.");
  ASSERT(checks[Interpolant::check_b], "Right boundary of the interval is not set.");
  ASSERT(checks[Interpolant::check_size], "Step is not set.");
  ASSERT(! ((user_tol == 0) && automatic_step), "Tolerance for automatic interpolation is not set.");
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
    func_diffn(NULL)
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

void Interpolant::create_nodes()
{
  step = (bound_b_-bound_a_)/size_;       //n_nodes = size+1;
  n_nodes = size_ + 1;
  
  x_vec.resize(n_nodes);       //nodes
  f_vec.resize(n_nodes);       //function values in the nodes
  df_vec.resize(n_nodes);      //function derivates values in the nodes
    
  //filling the vector x and f
  DiffValue value;
  double temp_x = bound_a_;
  for(unsigned int i = 0; i < size_; i++)
  {
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
    
  //DBGMSG("number_of_nodes = %d\n", n_nodes);
}

int Interpolant::interpolate()
{   
  check_all();  //checks if all needed parameters were set
  unsigned int result;
  
  //switch off statistics due to error computation
  //is switched on automatically by reset_stat()
  use_statistics = false;

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
  
  if(automatic_step)
  {
    DBGMSG("Maximum size of interpolation table: %d\n",max_size);
    unsigned int k=0;
    while(x_vec.size()-1 < max_size/2) 
    {
      //(this->*interpolate_func)();    //POSSIBLE WAY TO USE MORE KINDS OF INTERPOLATION
      interpolate_p1();
      if(user_tol < error_)
      {
        size_ *= 2;              //double the size (i.e. halve the step)
        DBGMSG("Interpolating: %d   size: %d \t error: %f\n",k, size_, error_);
        result = 1;   //tolerance has not been satisfied
      }
      else 
      {
        DBGMSG("Size of the table set automaticaly to: %d after %d cycles.\n", size_,k);
        DBGMSG("Error of the interpolation is: %f\n", error_);
        result = 0;   //interpolation OK
        break;
      }
      k++;
    }
    if(x_vec.size()-1 > max_size/2) 
    { DBGMSG("User defined tolerance %f has not been satisfied with size of interpolation table %d.\n",user_tol,size_); }
  }
  else 
  {
    //(this->*interpolate_func)();    //POSSIBLE WAY TO USE MORE KINDS OF INTERPOLATION
    interpolate_p1();
    result = 0;   //interpolation OK
  }
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
    create_nodes();
    p1_vec.resize(n_nodes-1);    //linear coeficients
    p1d_vec.resize(n_nodes-1);    //linear coeficients
    
    for(unsigned int i = 0; i < p1_vec.size(); i++)
    {
      p1_vec[i] = (f_vec[i+1] - f_vec[i]) / (x_vec[i+1] - x_vec[i]);
      p1d_vec[i] = (df_vec[i+1] - df_vec[i]) / (x_vec[i+1] - x_vec[i]);
    }
    
    FunctorBase<double>* norm = new NormW21(this);
    compute_error(norm);
    delete norm;
  
    //Writes the interpolation table.
    /*
    for(unsigned int i=0; i < p1_vec.size(); i++)
    {
      std::cout << "x: " << x_vec[i] << "\tf: " << f_vec[i] << "\tp1: " << p1_vec[i] << "\tdf: " << df_vec[i] << "\tp1d: " << p1d_vec[i] << std::endl;
    }
    //*/
  }
  
void Interpolant::compute_error(FunctorBase<double>* norm)
{
  double tot_err = 0; // total absolute error on <a,b>
  double p_err = 0;   // absolute error on x[i]-x[i+1]
  
  for(unsigned long i = 0; i < p1_vec.size(); i++ )
  { 
    p_err = std::sqrt(AdaptiveSimpson::AdaptSimpson(*norm,
                                                    x_vec[i], 
                                                    x_vec[i+1],
                                                    SIMPSON_TOLERANCE) );
    //DBGMSG("error on interval<%f,%f>: %f\n",x_vec[i],x_vec[i+1],p_err);
    tot_err += p_err;
  }  
  tot_err /= (bound_b_ - bound_a_);
    
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
  
  error_ = tot_err;       //returns total relative error
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


