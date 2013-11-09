#ifndef FUNCTORBASE_H
#define FUNCTORBASE_H

//include thirdparty/FADBAD++
#include "fadbad.h"
#include "badiff.h"
#include "tadiff.h"

using namespace fadbad;

///Class provides common functionality for functors.
/** This class provides methods for reading and setting functors
 * parameters.
 * 
 * User is adviced to create Enumeration to access parameters.
 * For example usage see class \ref FunctorBase.
 * 
 * One can also use \ref set_param_from_func
 * to copy all parameters from another functor.
 */
template<class Type>
class FunctorCommon
{
public:
  ///Constructor.
  FunctorCommon();

  ///Destructor.
  virtual ~FunctorCommon();
  
  ///Sets a functor's parameter.
  /**There is a vector for the functors parameters.
   * The vector of parameters is addressed by integer, but one can
   * use enum inside one's own functor to make it more convenient.
   * e.g. \verbatim set_param(MyFunctor::my_parameter, value) \endverbatim.
   * @param param_name is integer value (or enum value) to access member of the parameter vector
   * @param param_value is the new value of selected parameter
   */
  void set_param(unsigned int param_name, double param_value);
  
  ///Sets a functor's parameters from another functor.
  /** Copies all the parameters defined in the given functor.
   * @param func is pointer to a functor we want to copy from
   */
  template<class TType>
  void set_param_from_func(FunctorCommon<TType>* func);
  
  ///Returns parameter.
  /** @param param_name is integer value (or enum value) to access member of the parameter vector
   */
  double param(unsigned int param_name);
  
protected:
  std::vector<double> param_;
  
  ///Setting these friend classes to be able to access param_ vector.
  template<class U> friend class FunctorCommon;
};


///Abstract templated explicit functor class.
/** This class represents an explicit functor. Actual non-abstract functor
 * is obtained by deriving from this class and overriding the operator().
 * Parameters can be set and used by the methods of ancestor class FunctorCommon.
 * 
 * The template is used to enable computing derivates by FADBAD++ library. 
 * The virtual operator() needs to operate both with type double and FADBAD++ types, 
 * such as \verbatim B<double> \endverbatim and \verbatim T<double> \endverbatim.
 * 
 * Simple usage example for function \[x^3\] with three parameters: \code
template<class Type=double>
class Cubic : public FunctorBase<Type>
{
public:
  typedef enum{ p1, p2, p3
  } Parameters;

  Cubic(){}
  
  Cubic(double pp1, double pp2, double pp3)
  {
    this->set_param(p1, pp1);
    this->set_param(p2, pp2);
    this->set_param(p3, pp3);
  }
  
  virtual Type operator()(Type x)
  {
    return x*x*x * this->param(p1) + this->param(p2)/this->param(p3);
  }    
};
 * \endcode 
 */
template<class Type>
class FunctorBase : public FunctorCommon<Type>
{
public:
  ///Constructor.
  FunctorBase();

  ///Destructor.
  virtual ~FunctorBase();
  
  ///Virtual operator () with type @p Type.
  virtual Type operator()(Type x) = 0;
};



///Abstract templated implicit functor class.
/**This class represents an implicit functor.
 * It is used the same way as the class \ref FunctorBase,
 * only the operator() has two arguments.
 */
template<class Type>
class IFunctorBase : public FunctorCommon<Type>
{
public:
  ///Constructor.
  IFunctorBase();

  ///Destructor.
  virtual ~IFunctorBase();
  
  ///Virtual operator () with type @p Type.
  virtual Type operator()(Type x, Type y) = 0;
};

#endif  //FUNCTORBASE_H