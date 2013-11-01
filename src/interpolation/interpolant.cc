#include "functors_impl.hh"
#include "interpolant_impl.hh"
#include "adaptivesimpson.hh"

#include "system/xio.h"
//#include "system/sys_profiler.hh"

/********************************** InterpolantBase ******************************/
InterpolantBase::InterpolantBase()
  : degree(0),
    max_size(MAX_SIZE),
    automatic_step(false),
    error_(-1.0),
    extrapolation(InterpolantBase::functor),
    n_checks(4)
{
  reset_stat();
  checks.resize(n_checks);
}
  
InterpolantBase::~InterpolantBase() {}
  

double InterpolantBase::error()
{
  return error_;
}

void InterpolantBase::set_interval(const double& bound_a, const double& bound_b)
{
  ASSERT(bound_a!=bound_b,"Bounds overlap.");
  ASSERT(bound_a<bound_b,"a must be lower and b must be upper bound.");
  this->bound_a = bound_a;
  checks[check_a] = true;
  this->bound_b = bound_b;
  checks[check_b] = true;
}
  
void InterpolantBase::set_size(const unsigned int& size)
{
  ASSERT(size >0, "Size of interpolation table must be positive number!.");
  this->size = size;
  checks[Interpolant::check_size] = true;
  
  automatic_step = false;
}

void InterpolantBase::set_size_automatic(const double& user_tol, const unsigned int& init_size, const unsigned int& max_size)
{
  ASSERT(init_size >0, "Maximal size of interpolation table must be positive number!.");
  ASSERT(max_size >= init_size, "Maximal size of interpolation table is smaller than initial size.");
  this->user_tol = user_tol;
  this->max_size = max_size;
  size = init_size;
  automatic_step = true;
}

void InterpolantBase::set_extrapolation(InterpolantBase::ExtrapolationType extrapolation)
{
  this->extrapolation = extrapolation;
}

void InterpolantBase::reset_stat()
{
  interval_miss_a = 0;
  interval_miss_b = 0;
  total_calls = 0;
  max_a = 0;
  max_b = 0;
}

void InterpolantBase::check_and_reinterpolate()
{
  double percetange = 0.3;      //percentage of misses that is allowed
  bool reinterpolate = false;   //should we remake interpolation?
  if(interval_miss_a/total_calls > percetange)
  {
    bound_a = max_a;
    reinterpolate = true;
  }
  if(interval_miss_b/total_calls > percetange)
  {
    bound_b = max_b;
    reinterpolate = true;
  }
  
  if(reinterpolate)
  {
    interpolate(degree);
    reset_stat();
  }
}

  
bool InterpolantBase::check_all()
{
  if (!checks[Interpolant::check_functor]) std::cout << "Functor is not set." << std::endl;
  if (!checks[Interpolant::check_a]) std::cout << "Left boundary of the interval is not set." << std::endl;
  if (!checks[Interpolant::check_b]) std::cout << "Right boundary of the interval is not set." << std::endl;
  if (!checks[Interpolant::check_size]) std::cout << "Step is not set." << std::endl;
    
  bool res=true;
  for(unsigned int i=0; i<checks.size(); i++)
    res = res && checks[i];
    
  return res;
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
  
Interpolant::~Interpolant() {}

  
double Interpolant::f_diffn(const double& i_x, const unsigned int& n)
{
  ASSERT(n <= DIFF_N,"Not allowed to obtain n-th derivate. Check/change 'DIFF_N' in header file.");
  T<double> x,f;

  x = i_x;
  x[1] = 1;
  f = func_diffn->operator()(x);
  f.eval(DIFF_N);
      
  for(int i = 0; i<DIFF_N; i++)    //goes through the field of Taylor's coeficients
  {                           //and divide them by the factorial to get derivates
    f[i] = f[i] * this->fact(i);
  }
    
  return f[n];                //returns n-th derivate
}

void Interpolant::create_nodes()
{
  step = (bound_b-bound_a)/size;       //n_nodes = size+1;
  n_nodes = size + 1;
  
  x_vec.resize(n_nodes);       //nodes
  f_vec.resize(n_nodes);       //function values in the nodes
  df_vec.resize(n_nodes);      //function derivates values in the nodes
    
  //filling the vector x and f
  der value;
  double temp_x = bound_a;
  for(unsigned int i = 0; i < size; i++)                
  {
    temp_x = bound_a + step*i;
    value = f_diff(temp_x);
    x_vec[i] = temp_x;
    f_vec[i] = value.f;
    df_vec[i] = value.dfdx;
  }  
    
  //finish the interval
  x_vec[size] = bound_b;
  value = f_diff(bound_b);
  f_vec[size] = value.f;
  df_vec[size] = value.dfdx;
    
  //DBGMSG("number_of_nodes = %d\n", n_nodes);
}

int Interpolant::interpolate(unsigned int degree)
{
  this->degree = degree;        
  ASSERT(check_all(), "Parameters check did not pass. Some of the parameters were not set.");  
  unsigned int result;
  
  //selecting interpolation
  void (Interpolant::*interpolate_func)(void);
  switch(degree)
  {
    case 0: interpolate_func = &Interpolant::interpolate_p0; break; 
    case 1: interpolate_func = &Interpolant::interpolate_p1; break; 
    default: ASSERT(degree > 1, "Higher order interpolation is not available at the moment.");
  }
  
  if(automatic_step)
  {
    DBGMSG("Maximum size of interpolation table: %d\n",max_size);
    unsigned int k=0;
    while(x_vec.size()-1 < max_size/2) 
    {
      (this->*interpolate_func)();
      if(user_tol < error_)
      {
        size *= 2;              //double the size (i.e. halve the step)
        DBGMSG("Interpolating: %d   size: %d \t error: %f\n",k, size, error_);
        result = 1;   //tolerance has not been satisfied
      }
      else 
      {
        DBGMSG("Size of the table set automaticaly to: %d after %d cycles.\n", size,k);
        DBGMSG("Error of the interpolation is: %f\n", error_);
        result = 0;   //interpolation OK
        break;
      }
      k++;
    }
    if(x_vec.size()-1 > max_size/2) 
    { DBGMSG("User defined tolerance %f has not been satisfied with size of interpolation table %d.\n",user_tol,size); }
  }
  else 
  {
    (this->*interpolate_func)();
    result = 0;   //interpolation OK
  }
  return result;
}

  void Interpolant::interpolate_p0()
  { 
    create_nodes();
    val_ = &Interpolant::val_p0;      //pointer to the evalutation function
    diff_ = &Interpolant::diff_p0;      //pointer to the evalutation function
    
    Functor<double>* norm = new NormW21(this);
    compute_error(norm);      //sets error_
    delete norm;
  }
  
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
    
    val_ = &Interpolant::val_p1;      //pointer to the evalutation function
    diff_ = &Interpolant::diff_p1;    //pointer to the evalutation function
    
    Functor<double>* norm = new NormW21(this);
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
  
void Interpolant::compute_error(Functor<double>* norm)
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
  tot_err /= (bound_b - bound_a);
    
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
    fix_(InterpolantImplicit::no_fix)
  {
  }
  
InterpolantImplicit::~InterpolantImplicit() {}

void InterpolantImplicit::fix_variable(InterpolantImplicit::fix_var fix, const double& val)
{
  fix_ = fix;
  fix_val = val;
  func_u = new InterpolantImplicit::FuncExplicit<double>(*func,fix_,fix_val);
}

double InterpolantImplicit::f_val(const double& u)
{
  return func_u->operator()(u);
}


void InterpolantImplicit::interpolate_p0()
{
  if(func_u != NULL) 
    delete func_u;
  
  func_u = new FuncExplicit<double>(*func,fix_,fix_val);
}


