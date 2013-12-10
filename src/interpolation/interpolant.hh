#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "functors.hh"

/** This pair contains both the function and the first derivative values and is used
 * to as return value for functions of interpolation objects.
 * The reason is that the FADBAD++ library computes both function value and derivative value
 * at the same time. So when we use FADBAD++ in interpolant classes to compute 
 * derivates, we also return both values and we stick to it in interpolation evaluation also.
 * We use the form \verbatim pair<function_value, derivative_value> \endverbatim. 
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
  
/** Enumerates possible norms in which the interpolation error can be computed.
   */
  struct ErrorNorm
  {
    typedef enum {l2,      ///< L2 norm
                  lp,      ///< Lp norm
                  w21,     ///< W21 norm
                  wp1,     ///< Wp1 norm
                  max      ///< Maximum norm
    } Type;
  };
  
/** Enumerates types of variable fixation in implicit interpolant.
 */
  struct IFixVariable
  {
    typedef enum { fix_x,  ///< x variable will be fixed with given value 
                   fix_y,  ///< y variable will be fixed with given value 
                   no_fix  ///< no variable is fixed (used when InterpolantImplicit is created) 
    } Type;
  };
  
///Base class for interpolation.
/** This class serves as the common interface for interpolation classes.
 * It provides functions for setting the common parameters of an interpolation
 * such as interval, size of the table and extrapolation.
 * 
 * @see Interpolant
 * @see InterpolantImplicit
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

  ///Sets the type of norm used for computing estimate of the error of the interpolation.
  /** @param norm_type is the type of norm in which the error of interpolation is computed
   * @param p is the exponent used in norms \f$ L_p \f$ and \f$ W_p^1 \f$
   */
  void set_norm(ErrorNorm::Type norm_type=ErrorNorm::max, double p=2);
  
  /// Sets automatic choice of the size of the table.
  /** When @p interpolate is called than the size is increased (by smaller interval dividing) 
   * until the given tolerance or the maximum size of the interpolation table is reached.
   * @param user_tol is the tolerance which should the interpolation meet
   * @param init_size is the initial choice of table size
   * @param max_size is maximal size of the table that user allows
   */
  void set_size_automatic(double user_tol, unsigned int init_size, unsigned int max_size=default_max_size);
  
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
   * Then it compares the results with given percentage and accordingly changes (only streches) the interval and reinterpolate.
   * @param percentage is the percentage of evaluations outside the interval (on one side)
   */
  void check_stats_and_reinterpolate(double percentage=0.3);
  //@}
  
protected:
  double bound_a_,      ///< Left interval boundary.
         bound_b_,      ///< Right interval boundary.
         step,          ///< Chosen interpolation step.
         a_div_step;    ///< bound_ divided by step - precomputed value for evaluation
         
  unsigned int size_,   ///< Number of dividing intervals.
               n_nodes; ///< Number of nodes in the interval \f$(a,b)\f$.
               
  double user_tol;      ///< User set tolerance which is used during automatic step choice.
  unsigned int max_size;///< Maximal size of the interpolation table.
  bool automatic_size;  ///< Is true if step/size should be chosen automatically.
  ErrorNorm::Type norm_type;   ///< Type of norm used to compute error of the interpolation.
  double p;             ///< Exponent used in norms \f$ L_p \f$ and \f$ W_p^1 \f$ when computing error.
  
  double error_;        ///< Error of the interpolation.
  
  /** Extrapolation type - 'what' is evaluated outside the interpolation interval.
   * Default value after interpolant construction is  InterpolantBase::functor.
   */
  Extrapolation::Type extrapolation;     
  
  
  const static unsigned int n_derivatives;      ///< Defines how many derivatives we allow to be returned from Taylor's coeficients.
  const static unsigned int default_max_size;   ///< Default maximal size of the interpolation table.
  const static double simpson_tolerance;        ///< Tolerance in Adaptive Simpson intergration. @see AdaptiveSimpson
  
  ///Recursive factorial function (used in Taylor row expansion in n-th derivative computation).
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
  //bool use_statistics;       ///< Is true if statistics is checked (it must be switched of during error computation)
  EvalStatistics stats;        ///< Structure which keeps evaluation statistics. See \ref InterpolantBase::eval_statistics.
  
};


/// The main class for interpolation of functors.
/** This class is a wrapper of a functor (defined by class \ref FunctorBase) which provides
 * computation of derivates and interpolation of the functor. We will describe the functionality
 * in the following paragrahphs.
 *
 * ## Evaluation of a functor and its derivatives
 *
 * Beside the interpolation, the class provides calling of the functor itself.
 * All functions evaluating directly the functor starts with '`f_`', e.g. function \ref f_val computes the functor value.
 * 
 * We use the FADBAD++ library ([website](http://www.fadbad.com)) to compute derivatives of the functor. First derivatives
 * can be obtained by function \ref f_diff which uses backward automatic differentiation. Higher order derivatives can be 
 * obtained via \ref f_diffn which uses Taylor expansion method.
 * 
 * ## Interpolation
 * 
 * Sofar we provide only a  piecewise linear interpolation. Both the function and its derivative can be
 * interpolated. User has to provide the functor object to the constructor or later to the function \ref set_functor.
 * Then one must set at least the interval \f$ [a,b]\f$ by \ref set_interval and the size of the interpolation table 
 * by \ref set_size before calling \ref interpolate. Otherwise assert will crash the program in debug build.
 * 
 * When user wants the size to be chosen automatically then \ref set_size_automatic must be called at first.
 * An initial guess of size of the interpolation table and a tolerance which is assumed to be the maximum relative difference
 * between the function and the interpolant (sum of value and derivative difference) must be provided. 
 * Futher user can set the maximum size of the interpolation (default 1000) table which will not be exceeded and the type of norm
 * in which the estimate of the error of interpolation will be computed (enumerated in \ref ErrorNorm, default maximal norm).
 * See one of the following paragraphs to learn how we deal with the error of the interpolation.
 * 
 * Outside the given interval the interpolation is extrapolated. User can choose among several types of extrapolation
 * which are enumerated in \ref Extrapolation. Constant and linear extrapolation use the first and last interpolating polynomials
 * to compute the value. Functor type means that the functor itself will be evaluated outside the interval. Type of extrapolation
 * can be changed by \ref set_extrapolation function. Calling functor is the default type of extrapolation.
 * 
 * 
 * ### Error of the interpolation
 * 
 * We are interested in the maximum difference between the function value and its interpolation.
 * When interpolating also the derivative, we need to compute the difference between derivative and
 * its interpolation too.
 * Therefore we would like to compute something like the L-infinity norm, relative to function value and its derivative.
 * Let \f$ f(x), f'(x) \f$ be the function and its derivative and \f$ g(x), g'(x)\f$ interpolation of the function and interpolation 
 * of the derivative on the interval \f$ I \f$. We suggest an approximation of the infinity norm: 
 * \f[ 
 * \|f-g\|_\infty = \sup_{x\in I} F(x) = \sup_{x\in I} \left( \frac{|f(x)-g(x)|}{|f(x)| + tol} + \frac{|f'(x)-g'(x)|}{|f'(x)| + tol} \right)
 * \approx \max_{x\in J} \left( \frac{|f(x)-g(x)|}{|f(x)| + tol} + \frac{|f'(x)-g'(x)|}{|f'(x)| + tol} \right)
 * \f]
 * where \f$ J \f$ is set of points that lies in the middle of two neighboring nodes. \f$ tol \f$ is a given tolerance for zero (avoids zero division).
 * When computing the size of the interpolation table automatically, we compare the user tolerance \f$ TOL \f$ given by \ref set_size_automatic funtion 
 * directly with the estimated value of the supremum norm.
 * 
 * This approach is fast and when computing the size automatically, we use the computed values (middle points, function and derivative values) in next
 * interpolation. On the other hand this does not work good around zero points of the function where the value is small and the fraction becomes large or does
 * not converge when dividing the intervals (\f$ x^2 \f$ is a good example).
 * 
 * Therefore we provide also error estimation with \f$ L_p \f$ and \f$ W_p^1 \f$ norm. Then the user given tolerance is compared
 * as it is stated in the following inequation:
 * \f[
 * \left( \frac{\int \limits_I F^p(x) \,\rm{d} x}{|I|} \right)^{\frac{1}{p}} \leq TOL.
 * \f]
 * This approach is much more expensive but can give better results.
 * 
 * TODO: Improve error computation. Suggest more robust method. 
 * 
 * ### Evaluation statistics
 * 
 * We collect some statistics during the evaluation of the interpolant in the struct \ref InterpolantBase::EvalStatistics -- total number of calls, number of evaluations
 * inside the interpolation interval and number of evaluation outside the interval. One can use these statistics to decide whether it is necessary
 * to widen the interval and recompute the interpolation. Function \ref check_stats_and_reinterpolate is provided to do this automatically 
 * according to the given percentage of evaluation outside the interval.
 * 
 */

class Interpolant : public InterpolantBase
{
public:
  ///Default constructor.
  Interpolant();
  
  ///Constructor with functor setting.
  /** 
   * @param func is the pointer to functor
   * @param interpolate_derivative is true when derivate is also interpolated
   * @tparam Func is the functor type
   * @tparam Type is the template type of the functor (e.g. double)
   */
  template<template<class> class Func, class Type >
  Interpolant(Func<Type>* func, bool interpolate_derivative=false);
  
  ///Destructor.
  virtual ~Interpolant(void);
  
  ///@name Evaluation.
  //@{
  ///Returns interpolated value.
  /** @param x is the point at which we evaluate the interpolation
   */
  double val(double x);
  
  ///Do NOT use, only for testing purpose.
  /** TODO: After testing it can be removed and 
   * \ref val_p1 can be made private again.
   */
  double val_test(double x);
  
  ///Do NOT use, unless you are 100% sure. 
  /** This function evaluates the P1 interpolant at x.
   * It does not check the interval, does not provide extrapolation 
   * and does not collect statistics.
   * 
   * Can be used to POSSIBLY speed the evaluation just a little bit,
   * if you are absolutely sure that you evaluate the interpolant
   * only on the given interval and do not want to collect statistics.
   * 
   * Same can be done with \ref diff_p1 if it is made public.
   * 
   * Used in unit_test benchmark to compare with val function.
   */
  double val_p1(double x);

  ///Returns interpolated value of the derivation.
  /** @param x is the point at which we evaluate the interpolation
   */
  DiffValue diff(double x);
  
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
   * @param interpolate_derivative is true when derivate is also interpolated
   * @tparam Func is the functor type
   * @tparam Type is the template type of the functor (e.g. double)
   */
  template<template<class> class Func, class Type >
  void set_functor(Func<Type>* func, bool interpolate_derivative=false);
  //@}
  
protected:
  class FuncError_lp;
  class FuncError_wp1;
  
  FunctorBase<double>* func;                ///< Pointer to original functor with double type.
  FunctorBase<B<double> >* func_diff;       ///< Pointer to original functor with FADBAD B type.
  FunctorBase<T<double> >* func_diffn;      ///< Pointer to original functor with FADBAD T type.
  
  bool interpolate_derivative;  ///< Is true if we want to interpolate the derivative too.
  
  //std::vector<double> x_vec;    ///< Vector of nodes.
  std::vector<double> f_vec;    ///< Vector of function values at nodes.
  std::vector<double> df_vec;   ///< Vector of function derivatives values at nodes.
  //std::vector<double> p1_vec;   ///< Vector of linear coeficients of P1 interpolation.
  //std::vector<double> p1d_vec;  ///< Vector of linear coeficients of P1 interpolation.
  
  //Creates piecewise constant interpolation.
  //void interpolate_p0();
  
  ///Creates piecewise linear interpolation.
  //void interpolate_p1();
  
  ///Computes estimate of interpolation error in maximum norm.
  void compute_error(double tol, std::vector<double>& f, std::vector<double>& df);
  
  ///Computes estimate of interpolation error with given norm.
  void compute_error(double tol, double p, ErrorNorm::Type norm_type);
  
  void swap_middle_values(std::vector<double>& f, std::vector<double>& df);
  
  ///Creates table of nodes and function values.
  /** Creates vector of nodes according to the table size
   * and computes function and derivative values at the nodes.
   */
  void compute_values();

  ///Finds interval on which @p x lies.
  //unsigned int find_interval(double x);
  
  
  ///Function that evaluates the derivative of P1 interpolant at @p x.
  DiffValue diff_p1(double x);
};





class InterpolantImplicit : public InterpolantBase
{
public:
  
  ///@name Construction.
  //@{
  ///constructor
  InterpolantImplicit();
  
  ///Constructor with functor setting.
  /** 
   * @param func is the pointer to functor
   * @param interpolate_derivative is true when derivate is also interpolated
   * @tparam Func is the functor type
   * @tparam Type is the template type of the functor (e.g. double)
   */
  template<template<class> class Func, class Type >
  InterpolantImplicit(Func<Type>* func, bool interpolate_derivative=false);
  
  ///destructor
  virtual ~InterpolantImplicit(void);
  
  ///Sets the implicit functor.
  /** 
   * @param func is the pointer to implicit functor.
   * @tparam Func is the functor class.
   * @tparam Type is the template type of the functor (e.g. double)
   */
  template<template<class> class Func, class Type >
  void set_functor(Func<Type>* func, bool interpolate_derivative=false);
  
  //@}
  
  ///@name Evaluation.
  //@{
  
  /** Fixes the chosen variable and sets its fixed value.
   * @param fix is the chosen variable (no_fix, fix_x or fix_y)
   * @param value is the fixed value
   */
  void fix_variable(IFixVariable::Type fix, double value);
  
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
  virtual int interpolate();
  
  
  //@}
    
    
protected:
  
  template<class Type=double>
  class FuncExplicit;
  
  FuncExplicit<double>* func_u;
  IFunctorBase<double>* func;
  IFunctorBase<B<double> >* func_diff;
  IFunctorBase<T<double> >* func_diffn;
  
  bool interpolate_derivative;  ///< Is true if we want to interpolate the derivative too.
  
  Interpolant* explicit_interpolant;

  IFixVariable::Type fix_;
  double fix_val;
  
  ///@name Interpolation.
  //@{
  
  ///Creates piecewise linear interpolation.
  void interpolate_p1();
  
  //@}

};

#endif