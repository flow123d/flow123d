
#include "functors.hh"

#define SIMPSON_TOLERANCE 1e-10
#define MAX_SIZE 1000
#define MISS_PERCENTAGE 0.3      //percentage of interpolantion interval misses during evaluation that is allowed
#define STATISTIC_CHECK 1000      //frequance of checking the evaluation statistics


/** This pair contains both the function value and the first derivative.
 * We use it in interpolation to return both values. One of the reason
 * is that the FADBAD++ library computes both function value and derivative
 * at the same time.
 * We asume the form \verbatim pair<function_value, derivative> \endverbatim. 
 */
typedef std::pair<double,double> DiffValue;

///class Interpolant
/** Class can be templated by the functor (object with implemented operator ()), 
 * or there can still FuntorValueBase(), abstract class with virtual operator()
 * which would be passed to constructor.
 * It uses FADBAD 
 * library to obtain 1st derivative.
 */
class InterpolantBase
{
public:
  /** Enumerates possible kinds of interpolation.
   * @p functor means that the original functor will be called outside interpolation interval
   */
  typedef enum { constant, linear, functor 
  } ExtrapolationType;
  
  ///Structure that keeps statistics of evaluation.
  typedef struct {
  unsigned int interval_miss_a,         ///<< counts left misses of the interval
               interval_miss_b,         ///<< counts right misses of the interval
               total_calls;             ///<< counts total calls of evaluation
               
  double min,         ///<< minimal x for which the evaluation was called (outside the interval)
         max;         ///<< maximal x for which the evaluation was called (outside the interval)
  } eval_statistics;
  
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
  
  //Gettters
  eval_statistics statistics() const;
  double bound_a() const;
  double bound_b() const;
  unsigned int size() const;
  
  ///Sets the interpolation interval boundaries.
  void set_interval(double bound_a,double bound_b);
  
  ///Sets size of the interpolation table. It is also equal to the number of intervals.
  void set_size(unsigned int size);

  /** Sets automatic step choice.
   * When @p interpolate is called than divides the step until some tolerance or maximum
   * size of the interpolation table is reached.
   * @param user_tol is the tolerance which should the interpolation meet
   * @param init_size is the initial choice of table size
   * @param max_size is maximal size of the table that user allows
   */
  void set_size_automatic(double user_tol, unsigned int init_size, unsigned int max_size=MAX_SIZE);
  
  ///Sets the type of extrapolation. Functor type is default.
  void set_extrapolation(ExtrapolationType extrapolation);
  
  ///Creates piecewise interpolation with polynomials of selected degree.
  /** @return 0 if interpolation created; 1 if the tolerance has not been satisfied (in case of automatic choice of size)
   */
  virtual int interpolate() = 0;
  
  void check_and_reinterpolate();
  //@}
    
protected:
  ///Defines how many derivatives we allow to be returned from Taylor's coeficients.
  const static unsigned int n_derivatives;
  
  ///@name Interpolation.
  //@{
  double bound_a_,       ///<< Left interval boundary.
         bound_b_,       ///<< Right interval boundary.
         step;          ///<< Chosen interpolation step.
         
  unsigned int size_,    ///<< Number of dividing intervals.
               n_nodes; ///<< Number of nodes in the interval \[(a,b)\].
               
  double user_tol;      ///<< User set tolerance which is used during automatic step choice.
  unsigned int max_size; ///<< Maximal size of the interpolation table.
  bool automatic_step;  ///<< Is true if step/size should be chosen automatically.
  
  double error_;        ///<< Error of the interpolation. (Norm of difference divided by the lenght of interval.)
  
  ExtrapolationType extrapolation;      ///Extrapolation type.
  //@}
  
  ///@name Check.
  //@{
  ///Parameters setting check.
  enum {check_functor, check_a, check_b, check_size
  } check_type;

  ///Checks that the parameters are set before interpolation.
  std::vector<bool> checks;
  void check_all();
  //@}
  
  ///@name Evaluation statistics.
  //@{
  bool use_statistics;          ///<< true if statistics is checked (it must be switched of during error computation)
  eval_statistics stats;        ///<< structure which keeps evaluation statistics
  void reset_stat();            ///<< resets all measured statistics
  //@}
  
  ///Factorial (used in Taylor row expansion in n-th derivative computation @p diffn)
  long fact (long x);
};




class Interpolant : public InterpolantBase
{
public:
  ///@name Construction.
  //@{
  ///constructor
  Interpolant();
  
  ///Constructor with functor setting.
  /** 
   * @param func is the pointer to functor.
   * @tparam Func is the functor type
   * @tparam Type is the template type of the functor (e.g. double)
   */
  template<template<class> class Func, class Type >
  Interpolant(Func<Type>* func);
  
  ///destructor
  virtual ~Interpolant(void);
  
  
  ///Sets the functor.
  /** 
   * @param func is the pointer to functor.
   * @tparam Func is the functor type
   * @tparam Type is the template type of the functor (e.g. double)
   */
  template<template<class> class Func, class Type >
  void set_functor(Func<Type>* func);
  //@}
  
  
  ///@name Evaluation.
  //@{
  ///Returns interpolated value.
  /** @param x is the point at which we evaluate the interpolation
   */
  double val(double x);
  
  ///Returns interpolated value of the derivation.
  /** @param x is the point at which we evaluate the interpolation
   */
  DiffValue diff(double x);
  
  /** Returns interpolated n-th derivative.
    * @param x is the point at which we evaluate the original functor
    * @param n is the order of the derivative we want
    */
  //double diffn(double x, unsigned int n);
  
  ///Returns value of the original functor.
  /** @param x is the point at which we evaluate the original functor
   */
  double f_val(double x);
  
  ///Returns 1st derivative of original functor using FADBAD.
  /** @param x is the point at which we evaluate the original functor
   */
  DiffValue f_diff (double x);
  
  /** Returns n-th derivative of original functor using FADBAD.
    * Uses coeficients in Taylor's row.
    * @param x is the point at which we evaluate the original functor
    * @param n is the order of the derivative we want
    */
  double f_diffn(double x, unsigned int n);
  //@}
  

  ///@name Interpolation.
  //@{
  ///Creates piecewise interpolation with polynomials of selected degree.
  virtual int interpolate();
  //@}
    
    
protected:
  // Not actual norm.
  // Returns only power of difference of functor and its interpolant at point.
  class NormL2;
  class NormW21;
  
  FunctorBase<double>* func;                ///<< Pointer to original functor with double type.
  FunctorBase<B<double> >* func_diff;       ///<< Pointer to original functor with FADBAD type.
  FunctorBase<T<double> >* func_diffn;      ///<< Pointer to original functor with FADBAD type.
  
  ///@name Interpolation.
  //@{
  std::vector<double> x_vec;    ///<< Vector of nodes.
  std::vector<double> f_vec;    ///<< Vector of function values at nodes.
  std::vector<double> df_vec;   ///<< Vector of function derivatives values at nodes.
  std::vector<double> p1_vec;   ///<< Vector of linear coeficients of P1 interpolation.
  std::vector<double> p1d_vec;  ///<< Vector of linear coeficients of P1 interpolation.
  
  ///Creates piecewise constant interpolation.
  //void interpolate_p0();
  ///Creates piecewise linear interpolation.
  void interpolate_p1();
  
  ///Computes interpolation error with given norm.
  void compute_error(FunctorBase<double>* norm);
  
  /** Creates vector of nodes according to the table size
   * and computes function values and derivatives at the nodes.
   */
  void create_nodes();

  ///Finds interval on which @p x lies.
  unsigned int find_interval(double x);

  ///Function that evaluates the P1 interpolant at @p x.
  double val_p1(double x);

  ///Function that evaluates the derivative of P1 interpolant at @p x.
  DiffValue diff_p1(double x);

  /* CONSTANT INTERPOLATION
  double val_p0(double x);
  der diff_p0(double x);
  */
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
  
  ///Sets the implicit functor.
  /** 
   * @param func is the pointer to implicit functor.
   * @tparam Func is the functor class.
   * @tparam Type is the template type of the functor (e.g. double)
   */
  template<template<class> class Func, class Type >
  void set_functor(Func<Type>* func);
  
  //@}
  
  ///@name Evaluation.
  //@{
  
  /** Fixes the chosen variable and sets its fixed value.
   * @param fix is the chosen variable (no_fix, fix_x or fix_y)
   * @param value is the fixed value
   */
  void fix_variable(InterpolantImplicit::fix_var fix, const double& val);
  
  ///Returns interpolated value.
  /** @param u is the point at which we evaluate the interpolation
   */
  double val(double u);
  
  ///Returns interpolated value of the derivation.
  /** @param u is the point at which we evaluate the interpolation
   */
  DiffValue diff(double u);
  
  ///Returns interpolated value of the derivation.
  /** @param x is the point at which we evaluate the interpolation
   */
  //double diffn(const double& x, const unsigned int& n);
  
  ///Returns value of the original functor.
  /** @param x is function variable.
   *  @param y is function variable.
   */
  double f_val(double x, double y);
  
  ///Returns value of the original functor when one of the variables is fixed.
  /** @param u is function variable (the one not fixed).
   */
  double f_val(double u);
  
  ///Returns 1st derivative of original functor using FADBAD.
  /** @param x is function variable.
   *  @param y is function variable.
   */
  DiffValue f_diff (double x, double y);
  
  /** Returns n-th derivative of original functor using FADBAD.
    * Uses coeficients in Taylor's row.
    * @param u is the point at which we evaluate the original functor
    * @param n is the order of the derivative we want
    */
  double f_diffn(double u, unsigned int n);
  //@}
  

  ///@name Interpolation.
  //@{
  
  ///Creates piecewise interpolation with polynomials of selected degree.
  virtual int interpolate() {return 5;}
  
  ///Creates piecewise linear interpolation.
  void interpolate_p1();
  //@}
    
    
protected:
  
  template<class Type=double>
  class FuncExplicit;
  
  FuncExplicit<double>* func_u;
  IFunctorBase<double>* func;
  IFunctorBase<B<double> >* func_diff;
  IFunctorBase<T<double> >* func_diffn;
  
  Interpolant* explicit_interpolant;

  fix_var fix_;
  double fix_val;
  
  ///@name Interpolation.
  //@{
  
  //@}

};

