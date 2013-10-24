
#include "functors.hh"

///class Interpolant
/** Class can be templated by the functor (object with implemented operator ()), 
 * or there can still FuntorValueBase(), abstract class with virtual operator()
 * which would be passed to constructor.
 * It uses FADBAD 
 * library to obtain 1st derivate.
 */

class Interpolant
{
public:
  ///@name Construction.
  //@{
  ///constructor
  Interpolant();
  
  ///destructor
  virtual ~Interpolant(void);
  
  ///Sets the functor.
  /** 
   * @param func is the pointer to functor.
   * @tparam Func is the functor type
   * @tparam Type is the template type of the functor (e.g. double)
   */
  template<template<class> class Func, class Type >
  void set_functor(Functor<Type>* func);
  
  ///Specialization of the setter for Type=double
  template<template<class> class Func>
  void set_functor(Functor<double>* func);
  //@}
  
  ///@name Evaluation.
  //@{
  ///Returns interpolated value.
  /** @param i_x is the point at which we evaluate the interpolation
   */
  double val(const double& i_x);
  
  ///Returns interpolated value of the derivation.
  /** @param i_x is the point at which we evaluate the interpolation
   */
  der diff(const double& i_x);
  
  ///Returns interpolated value of the derivation.
  /** @param i_x is the point at which we evaluate the interpolation
   */
  //double diffn(const double& i_x, const unsigned int& n);
  
  ///Returns value of the original functor.
  /** @param i_x is the point at which we evaluate the original functor
   */
  double f_val(const double& i_x);
  
  ///Returns 1st derivate of original functor using FADBAD.
  /** @param i_x is the point at which we evaluate the original functor
   */
  der f_diff ( const double& i_x );
  
  /** Returns n-th derivate of original functor using FADBAD.
    * Uses coeficients in Taylor's row.
    * @param i_x is the point at which we evaluate the original functor
    * @param n is the order of the derivate we want
    */
  double f_diffn(const double& i_x, const unsigned int& n);
  //@}
  

  ///@name Interpolation.
  //@{
    ///Sets the interpolation interval.
    void set_interval(const double&a, const double& b);
    
    void set_size(const unsigned int& size);
    double set_size_automatic(const double& tol,const unsigned int& max_size=0);
    
    ///Creates piecewise constant interpolation.
    void interpolate_p0(double& interpolation_error);
    
    ///Creates piecewise linear interpolation.
    void interpolate_p1(double& interpolation_error);
  //@}
    
    
protected:
  Functor<double>* func;
  Functor<B<double> >* func_diff;
  Functor<T<double> >* func_diffn;
  
  ///@name Interpolation.
  //@{
  double bound_a,       ///<< Left interval boundary.
         bound_b,       ///<< Right interval boundary.
         step;          ///<< Chosen interpolation step.
         
  unsigned int size,    ///<< Number of dividing intervals.
               n_nodes; ///<< Number of nodes in the interval \[(a,b)\].
  
  std::vector<double> x_vec;    ///<< Vector of nodes.
  std::vector<double> f_vec;    ///<< Vector of function values at nodes.
  std::vector<double> df_vec;   ///<< Vector of function derivates values at nodes.
  std::vector<double> p1_vec;   ///<< Vector of linear coeficients of P1 interpolation.
  std::vector<double> p1d_vec;   ///<< Vector of linear coeficients of P1 interpolation.
  
  void create_nodes();
  double compute_error();
  
  unsigned int find_interval(const double& i_x);
  double (Interpolant::*val_)(const double&);
  der (Interpolant::*diff_)(const double&);
  double val_p0(const double& i_x);
  double val_p1(const double& i_x);
  der diff_p0(const double& i_x);
  der diff_p1(const double& i_x);
  //@}
  
  ///@name Check.
  //@{
  ///Parameters setting check.
  enum {check_functor, check_a, check_b, check_size
  } check_type;
  ///Number of checks.
  const double n_checks;
  ///Checks that the parameters are set before interpolation.
  std::vector<bool> checks;
  bool check_all();
  //@}
  
  ///@name Evaluation statistics.
  //@{
  ///Parameters setting check.
  unsigned int interval_miss_a,
               interval_miss_b,
               interval_hits;
               
  //@}
  
  ///Factorial (used in Taylor row expansion in n-th derivate computation @p diffn)
  long fact ( long x );
};




/******************************** Interpolant template implementation ********************************/
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