#ifndef FUNCTORBASE_H
#define FUNCTORBASE_H

//include thirdparty/FADBAD++
#include "fadbad.h"
#include "badiff.h"
#include "tadiff.h"

using namespace fadbad;

///Abstract templated class with virtual operator ()
template<class Type>
class FunctorBase
{
public:
  ///Constructor.
  FunctorBase();

  ///Destructor.
  virtual ~FunctorBase();
  
  ///Sets a functor's parameter.
  /**There is a vector for data needed inside functor (parameters).
   * The vector of parameters is addressed by integer, but one can
   * use enum inside one's own functor to make it more convenient.
   * e.g. \verbatim set_param(MyFunctor::my_parameter, value) \endverbatim.
   */
  void set_param(const unsigned int& param_name, const double& param_value);
  
  ///Sets a functor's parameters from another functor.
  /** Copies all the parameters defined in the given functor.
   */
  template<class TType>
  void set_param_from_func(FunctorBase<TType>* func);
  
  ///Returns parameter.
  double param(const unsigned int& param_name);
  
  /** Returns size of the vector for the parameters.
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

  ///Destructor.
  virtual ~FunctorImplicit();
  
  ///Virtual operator () with type @p Type.
  virtual Type operator()(const Type& x, const Type& y) = 0;
};

#endif  //FUNCTORBASE_H