#ifndef FUNCTORBASE_H
#define FUNCTORBASE_H

//include thirdparty/FADBAD++
#include "fadbad.h"
#include "badiff.h"
#include "tadiff.h"


///Defines how many derivates can be returned from Taylor's coeficients.
#define DIFF_N 10
//#define SIMPSON_TOLERANCE 1e-9

using namespace fadbad;

/// Structure der
/** Structure with derivate and value, that is return by Diff()
  */
struct der
{
  double f;
  double dfdx;
};

///Abstract templated class with virtual operator ()
template<class Type>
class FunctorBase
{
public:
  ///Constructor.
  FunctorBase();
  
  ///Copy constructor - templated.
  /** Enables doing copies like: \verbatim FunctorBase<B<double> >(FunctorBase<double>)\endverbatim
   * which is used in Interpolant's method @p Interpolant::set_functor.
   */
  template<class TType>
  FunctorBase(FunctorBase<TType>& func);

  ///Destructor.
  virtual ~FunctorBase();
  
  ///Sets a functor's parameter.
  /**There is a vector for data needed inside functor (parameters).
   * The vector of parameters is addressed by integer, but one can
   * use enum inside one's own functor to make it more convenient.
   * e.g. \verbatim set_param(MyFunctor::my_parameter, value) \endverbatim.
   */
  void set_param(const unsigned int& param_name, const double& param_value);
  
  ///Returns parameter.
  double param(const unsigned int& param_name);
  
  /** Returns size of the vector for the parameters.
   * Not actual number of used parameters! 
   */
  unsigned int n_param();
  
protected:
  std::vector<double> param_;
};


///Abstract templated class with virtual operator ()
template<class Type>
class Functor : public FunctorBase<Type>
{
public:
  ///Constructor.
  Functor();
  
  ///Copy constructor - templated.
  /** Enables doing copies like: \verbatim FunctorBase<B<double> >(FunctorBase<double>)\endverbatim
   * which is used in Interpolant's method @p Interpolant::set_functor.
   */
  template<class TType>
  Functor(FunctorBase<TType>& func);

  ///Destructor.
  virtual ~Functor();
  
  ///Virtual operator () with type @p Type.
  virtual Type operator()(const Type& x) = 0;
};



///Abstract templated class with virtual operator ()
template<class Type>
class FunctorImplicit : public FunctorBase<Type>
{
public:
  ///Constructor.
  FunctorImplicit();
  
  ///Copy constructor - templated.
  /** Enables doing copies like: \verbatim FunctorBase<B<double> >(FunctorBase<double>)\endverbatim
   * which is used in Interpolant's method @p Interpolant::set_functor.
   */
  template<class TType>
  FunctorImplicit(FunctorImplicit<TType>& func);

  ///Destructor.
  virtual ~FunctorImplicit();
  
  ///Virtual operator () with type @p Type.
  virtual Type operator()(const Type& x, const Type& y) = 0;
};

#endif  //FUNCTORBASE_H