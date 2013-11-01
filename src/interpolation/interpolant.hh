
#include "functors.hh"

#define SIMPSON_TOLERANCE 1e-10
#define MAX_SIZE 1000

///class Interpolant
/** Class can be templated by the functor (object with implemented operator ()), 
 * or there can still FuntorValueBase(), abstract class with virtual operator()
 * which would be passed to constructor.
 * It uses FADBAD 
 * library to obtain 1st derivate.
 */
class InterpolantBase
{
public:
  typedef enum { constant, linear, functor 
  } ExtrapolationType;
  
  ///@name Construction.
  //@{
  ///constructor
  InterpolantBase();
  
  ///destructor
  virtual ~InterpolantBase(void);
 
  ///@name Interpolation.
  //@{
  ///Returns error of the interpolation.
  ///It is equal to the chosen norm of the difference divided by the lenght of interval.
  ///Returns -1.0 if the interpolation has not been computed yet.
  double error();
  
  ///Sets the interpolation interval.
  void set_interval(const double&bound_a, const double& bound_b);
  
  void set_size(const unsigned int& size);
  void set_size_automatic(const double& user_tol,const unsigned int& init_size, const unsigned int& max_size=MAX_SIZE);
  
  void set_extrapolation(ExtrapolationType extrapolation);
  
  ///Creates piecewise interpolation with polynomials of selected degree.
  virtual int interpolate(unsigned int degree) = 0;
  //@}
    
protected:
  ///@name Interpolation.
  //@{
  double bound_a,       ///<< Left interval boundary.
         bound_b,       ///<< Right interval boundary.
         step;          ///<< Chosen interpolation step.
         
  unsigned int size,    ///<< Number of dividing intervals.
               n_nodes, ///<< Number of nodes in the interval \[(a,b)\].
               degree;
               
  double user_tol;
  unsigned int max_size;
  bool automatic_step;
  
  double error_;        ///<< Error of the interpolation. (Norm of difference divided by the lenght of interval.)
  
  ExtrapolationType extrapolation;
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
               total_calls;   
               
  double max_a,
         max_b;
         
  void reset_stat();
  void check_and_reinterpolate();
  //@}
  
  ///Factorial (used in Taylor row expansion in n-th derivate computation @p diffn)
  long fact ( long x );
};




class Interpolant : public InterpolantBase
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
  ///Creates piecewise interpolation with polynomials of selected degree.
  virtual int interpolate(unsigned int degree);
  //@}
    
    
protected:
  class NormL2;
  class NormW21;
  
  Functor<double>* func;
  Functor<B<double> >* func_diff;
  Functor<T<double> >* func_diffn;
  
  ///@name Interpolation.
  //@{
  std::vector<double> x_vec;    ///<< Vector of nodes.
  std::vector<double> f_vec;    ///<< Vector of function values at nodes.
  std::vector<double> df_vec;   ///<< Vector of function derivates values at nodes.
  std::vector<double> p1_vec;   ///<< Vector of linear coeficients of P1 interpolation.
  std::vector<double> p1d_vec;  ///<< Vector of linear coeficients of P1 interpolation.
  
  ///Creates piecewise constant interpolation.
  void interpolate_p0();
  ///Creates piecewise linear interpolation.
  void interpolate_p1();
  
  void compute_error(Functor<double>* norm);
  
  void create_nodes();
  unsigned int find_interval(const double& i_x);
  double (Interpolant::*val_)(const double&);
  der (Interpolant::*diff_)(const double&);
  double val_p0(const double& i_x);
  double val_p1(const double& i_x);
  der diff_p0(const double& i_x);
  der diff_p1(const double& i_x);
  //@}
};






class InterpolantImplicit : public InterpolantBase
{
public:
  typedef enum { fix_x, fix_y, no_fix
  } fix_var;
  
  ///@name Construction.
  //@{
  ///constructor
  InterpolantImplicit();
  
  ///destructor
  virtual ~InterpolantImplicit(void);
  
  ///Sets the functor.
  /** 
   * @param func is the pointer to functor.
   * @tparam Func is the functor type
   * @tparam Type is the template type of the functor (e.g. double)
   */
  template<template<class> class Func, class Type >
  void set_functor(FunctorImplicit<Type>* func);
  
  ///Specialization of the setter for Type=double
  template<template<class> class Func>
  void set_functor(FunctorImplicit<double>* func);
  //@}
  
  ///@name Evaluation.
  //@{
  
  void fix_variable(InterpolantImplicit::fix_var fix, const double& val);
  
  ///Returns interpolated value.
  /** @param u is the point at which we evaluate the interpolation
   */
  double val(const double& u);
  
  ///Returns interpolated value of the derivation.
  /** @param u is the point at which we evaluate the interpolation
   */
  der diff(const double& u);
  
  ///Returns interpolated value of the derivation.
  /** @param x is the point at which we evaluate the interpolation
   */
  //double diffn(const double& x, const unsigned int& n);
  
  ///Returns value of the original functor.
  /** @param x is function variable.
   *  @param y is function variable.
   */
  double f_val(const double& x, const double& y);
  
  ///Returns value of the original functor when one of the variables is fixed.
  /** @param u is function variable (the one not fixed).
   */
  double f_val(const double& u);
  
  ///Returns 1st derivate of original functor using FADBAD.
  /** @param x is function variable.
   *  @param y is function variable.
   */
  der f_diff ( const double& x, const double& y);
  
  /** Returns n-th derivate of original functor using FADBAD.
    * Uses coeficients in Taylor's row.
    * @param u is the point at which we evaluate the original functor
    * @param n is the order of the derivate we want
    */
  double f_diffn(const double& u, const unsigned int& n);
  //@}
  

  ///@name Interpolation.
  //@{
  
  ///Creates piecewise interpolation with polynomials of selected degree.
  int interpolate(unsigned int degree) {return 5;}
  
  ///Creates piecewise constant interpolation.
  void interpolate_p0();
  
  ///Creates piecewise linear interpolation.
  void interpolate_p1();
  //@}
    
    
protected:
  
  template<class Type=double>
  class FuncExplicit;
  
  FuncExplicit<double>* func_u;
  FunctorImplicit<double>* func;
  FunctorImplicit<B<double> >* func_diff;
  FunctorImplicit<T<double> >* func_diffn;
  
  fix_var fix_;
  double fix_val;
  
  ///@name Interpolation.
  //@{
  std::vector<double> x_vec;    ///<< Vector of nodes.
  std::vector<double> f_vec;    ///<< Vector of function values at nodes.
  std::vector<double> df_vec;   ///<< Vector of function derivates values at nodes.
  std::vector<double> p1_vec;   ///<< Vector of linear coeficients of P1 interpolation.
  std::vector<double> p1d_vec;  ///<< Vector of linear coeficients of P1 interpolation.
  
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

};

