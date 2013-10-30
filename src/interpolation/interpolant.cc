#include "functors_impl.hh"
#include "interpolant_impl.hh"
#include "adaptivesimpson.hh"

#include "system/xio.h"
//#include "system/sys_profiler.hh"

/********************************** InterpolantBase ******************************/
InterpolantBase::InterpolantBase()
  : n_checks(4),
    interval_miss_a(0),
    interval_miss_b(0),
    interval_hits(0)
{
  checks.resize(n_checks);
}
  
InterpolantBase::~InterpolantBase() {}
  

void InterpolantBase::set_interval(const double& a, const double& b)
{
  ASSERT(a!=b,"Bounds overlap.");
  ASSERT(a<b,"a must be lower and b must be upper bound.");
  this->bound_a = a;
  checks[check_a] = true;
  this->bound_b = b;
  checks[check_b] = true;
}
  
void InterpolantBase::set_size(const unsigned int& size)
{
  this->size = size;
  checks[Interpolant::check_size] = true;
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
    
  DBGMSG("number_of_nodes = %d\n", n_nodes);
}

  void Interpolant::interpolate_p0(double& interpolation_error)
  {
    ASSERT(check_all(), "Parameters check did not pass. Probably some of the parameters were not set.");
    create_nodes();
    val_ = &Interpolant::val_p0;      //pointer to the evalutation function
    diff_ = &Interpolant::diff_p0;      //pointer to the evalutation function
  }
  
  void Interpolant::interpolate_p1(double& interpolation_error)
  {
    ASSERT(check_all(), "Parameters check did not pass. Probably some of the parameters were not set.");
    create_nodes();
    p1_vec.resize(n_nodes-1);    //linear coeficients
    p1d_vec.resize(n_nodes-1);    //linear coeficients
    
    for(unsigned int i = 0; i < n_nodes; i++)
    {
      p1_vec[i] = (f_vec[i+1] - f_vec[i]) / (x_vec[i+1] - x_vec[i]);
      p1d_vec[i] = (df_vec[i+1] - df_vec[i]) / (x_vec[i+1] - x_vec[i]);
    }
    
    val_ = &Interpolant::val_p1;      //pointer to the evalutation function
    diff_ = &Interpolant::diff_p1;    //pointer to the evalutation function
    
    /*
    //Writes the interpolation table.
    for(unsigned int i=0; i<x_vec.size(); i++)
    {
      std::cout << "x: " << x_vec[i] << "\tf: " << f_vec[i] << "\tp1: " << p1_vec[i] << std::endl;
    }
    //*/
  }
  
  double Interpolant::compute_error()
  {
    double tot_err = 0; // absolute error on <a,b>
    /*
  ErrorNum p_err;       // absolute/relative polynomial error on its interval
      
  p_err.i = 0;
  p_err.err = 0;
        
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
  
  /* writes relative and absolute total error
  std::cout << "\trelative err = "
    << tot_err/(g->GetA()-g->GetB()) 
    << "\tabsolute err = " << tot_err << std::endl;
  //*/
  
    return tot_err;       //returns absolute total error
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
  func_u = new FuncExplicit<double>(*func,fix_,fix_val);
}

double InterpolantImplicit::f_val(const double& u)
{
  return func_u->operator()(u);
}


void InterpolantImplicit::interpolate_p0(double& interpolation_error)
{
  if(func_u != NULL) 
    delete func_u;
  
  func_u = new FuncExplicit<double>(*func,fix_,fix_val);
}


