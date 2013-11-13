
#include "functors.hh"

#define SIMPSON_TOLERANCE 1e-10
#define MAX_SIZE 1000
#define STATISTIC_CHECK 1000      //frequency of checking the evaluation statistics


/** This pair contains both the function and the first derivative values.
 * The reason is that the FADBAD++ library computes both function value and derivative value
 * at the same time. So when we use FADBAD++ in interpolant classes to compute 
 * derivates, we also return both values.
 * We assume the form \verbatim pair<function_value, derivative_value> \endverbatim. 
 */
typedef std::pair<double,double> DiffValue;


/** Enumerates possible kinds of extrapolation.
   */
  struct Extrapolation
  {
    typedef enum {constant,      ///< means that the values at the boundary points will be used outside the interval_miss_a
                  linear,        ///< means that the linear approximation on the first and last subintervals will be used outside the interval 
                  functor        ///< means that the original functor will be called outside interpolation interval
    } Type;
  };
  
  
///Base class for interpolation.
/** This class serves as the common interface for interpolation classes.
 * It provides classes for setting the basic parameters of an interpolation
 * such as interval, size of the table and extrapolation.
 */
class InterpolantBase
{
public:
  ///Default constructor.
  InterpolantBase();
  
  ///Destructor.
  virtual ~InterpolantBase();
  
  ///@name Interpolation.
  //@{
  ///Returns error of the interpolation.
  /**It is equal to the chosen norm of the difference divided by the lenght of interval.
   * Returns -1.0 if the interpolation has not been computed yet.
   */
  double error();
  
  //Gettters
  double bound_a() const;               ///< Returns left boundary of the interval.
  double bound_b() const;               ///< Returns right boundary of the interval.
  unsigned int size() const;            ///< Returns the size of the interpolation table.
  
  ///Sets the interpolation interval boundaries.
  /** @param bound_a is the left boundary of the interval
   *  @param bound_b is the right boundary of the interval
   */
  void set_interval(double bound_a,double bound_b);
  
  ///Sets size of the interpolation table. It is also equal to the number of intervals.
  /** @param size is the new size of the interpolation table
   */
  void set_size(unsigned int size);

  /// Sets automatic choice of the size of the table.
  /** When @p interpolate is called than the size is increased (by smaller interval dividing) 
   * until the given tolerance or the maximum size of the interpolation table is reached.
   * @param user_tol is the tolerance which should the interpolation meet
   * @param init_size is the initial choice of table size
   * @param max_size is maximal size of the table that user allows
   */
  void set_size_automatic(double user_tol, unsigned int init_size, unsigned int max_size=MAX_SIZE);
  
  ///Sets the type of extrapolation. Functor type is default.
  /** @param extrapolation is the type of extrapolation (defined by enumeration @p ExtrapolationType)
   */
  void set_extrapolation(Extrapolation::Type extrapolation);
  
  ///Creates piecewise polynomial interpolation.
  /** @return 
   * - 0 if interpolation has been created 
   * - 1 if the tolerance has not been satisfied (in case of automatic choice of size)
   */
  virtual int interpolate() = 0;
  //@}
  
  /** @name Evaluation statistics.
   * These members are used recording the statistics of the evaluation.
   */
  //@{
  ///Structure that keeps statistics of evaluation.
  typedef struct {
  unsigned int interval_miss_a,         ///< counts left misses of the interval
               interval_miss_b,         ///< counts right misses of the interval
               total_calls;             ///< counts total calls of evaluation
               
  double min,         ///< minimal x for which the evaluation was called outside the interval (initially is equal the left boundary)
         max;         ///< maximal x for which the evaluation was called outside the interval ((initially is equal the right boundary)
  } EvalStatistics;
  
  EvalStatistics statistics() const;   ///< Returns structure with evaluation statistics.
  
  ///Resets all measured statistics.
  void reset_stat();     
  
  ///Can be called to check automatically the evaluation statistics and possibly reinterpolate.
  /** The function computes ratio of evaluations outside the interpolation interval - on each side of the interval. 
   * Then it compares the results with given percentage and accordingly changes (streches) the interval a reinterpolate.
   * @param percentage is the percentage of evaluations outside the interval (on one side)
   */
  void check_stats_and_reinterpolate(double percentage=0.3);
  //@}
  
protected:
  double bound_a_,      ///< Left interval boundary.
         bound_b_,      ///< Right interval boundary.
         step;          ///< Chosen interpolation step.
         
  unsigned int size_,   ///< Number of dividing intervals.
               n_nodes; ///< Number of nodes in the interval \f$(a,b)\f$.
               
  double user_tol;      ///< User set tolerance which is used during automatic step choice.
  unsigned int max_size;///< Maximal size of the interpolation table.
  bool automatic_step;  ///< Is true if step/size should be chosen automatically.
  
  double error_;        ///< Error of the interpolation. (Norm of difference divided by the lenght of interval.)
  
  /** Extrapolation type - 'what' is evaluated outside the interpolation interval.
   * Default value after interpolant construction is  InterpolantBase::functor.
   */
  Extrapolation::Type extrapolation;     
  
  ///Defines how many derivatives we allow to be returned from Taylor's coeficients.
  const static unsigned int n_derivatives;
  ///Factorial (used in Taylor row expansion in n-th derivative computation).
  long fact (long x);

  
  ///@name Internal check system.
  /// This is internal check system which checks the values of parameters before making interpolation.
  //@{
  struct Check{
    ///Enumerates parameters that must be set before creation of the interpolant.
    enum {functor, bound_a, bound_b, size
    } Type;
  };
  ///Vector of boolean values telling us which parameters are set or not.
  std::vector<bool> checks;
  ///Checks that the parameters are set before interpolation.
  void check_all();
  //@}
  
  
  //currently not used
  //bool use_statistics;        ///< Is true if statistics is checked (it must be switched of during error computation)
  EvalStatistics stats;        ///< Structure which keeps evaluation statistics. See InterpolantBase::eval_statistics.
  
};



class Interpolant : public InterpolantBase
{
public:
  ///Default constructor.
  Interpolant();
  
  ///Constructor with functor setting.
  /** 
   * @param func is the pointer to functor
   * @tparam Func is the functor type
   * @tparam Type is the template type of the functor (e.g. double)
   */
  template<template<class> class Func, class Type >
  Interpolant(Func<Type>* func);
  
  ///Destructor.
  virtual ~Interpolant(void);
  
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
  
  ///Sets the functor.
  /** Can be used when the functor is not set in the constructor.
   * @param func is the pointer to functor
   * @tparam Func is the functor type
   * @tparam Type is the template type of the functor (e.g. double)
   */
  template<template<class> class Func, class Type >
  void set_functor(Func<Type>* func);
  //@}
    
    
protected:
  class NormL2;
  class NormW21;
  class NormWp1;
  class NormFunc;
  
  FunctorBase<double>* func;                ///< Pointer to original functor with double type.
  FunctorBase<B<double> >* func_diff;       ///< Pointer to original functor with FADBAD B type.
  FunctorBase<T<double> >* func_diffn;      ///< Pointer to original functor with FADBAD T type.
  
  std::vector<double> x_vec;    ///< Vector of nodes.
  std::vector<double> f_vec;    ///< Vector of function values at nodes.
  std::vector<double> df_vec;   ///< Vector of function derivatives values at nodes.
  std::vector<double> p1_vec;   ///< Vector of linear coeficients of P1 interpolation.
  std::vector<double> p1d_vec;  ///< Vector of linear coeficients of P1 interpolation.
  
  ///Creates piecewise constant interpolation.
  //void interpolate_p0();
  ///Creates piecewise linear interpolation.
  void interpolate_p1();
  
  double extrapolate_val(bool interval, double x);
  DiffValue extrapolate_diff(bool interval, double x);
  
  ///Computes interpolation error with given norm.
  void compute_error(FunctorBase<double>* norm);
  
  ///Creates table of nodes and function values.
  /** Creates vector of nodes according to the table size
   * and computes function and derivative values at the nodes.
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

