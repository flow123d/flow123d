/*
 * time_governor_test.cpp
 *
 *  Created on: May 20, 2011
 *      Author: jb
 */

#include <gtest/gtest.h>

#include "interpolation/functors.hh"
#include "interpolation/interpolant_impl.hh"

#define EQUAL(a,b) INPUT_CHECK( (a) == (b), #a": %f and "#b":%f differs\n",a,b);


//Example functor f(x)=x^3 
template<class Type=double>
class MyFunction_x3 : public Functor<Type>
{
public:
    
  typedef enum{ p1, p2, p3
  } Parameters;

  ///Constructor.
  MyFunction_x3(){}
  
  template<class TType>
  MyFunction_x3(Functor<TType>& func) : Functor<Type>(func){};
  
  virtual Type operator()(const Type& x)
  {
    return x*x*x;
  }    
};

/**
 * Test for Functors and interpolant creation only
 */
TEST (Functors, functors)
{
  MyFunction_x3<double> my_func;
  MyFunction_x3<> my_func2;
  
  my_func.set_param(MyFunction_x3<double>::p1,1.0);
  my_func.set_param(MyFunction_x3<>::p2,2.0);
  my_func.set_param(MyFunction_x3<>::p3,3.0);
  
  Interpolant* interpolant = new Interpolant();
  interpolant->set_functor<MyFunction_x3, double>(&my_func);
  
  Interpolant* interpolant2 = new Interpolant();
  interpolant2->set_functor<MyFunction_x3>(&my_func2);
  
  //params
  EQUAL(my_func.get_param(2), 3.0);
  EQUAL(my_func.get_param(MyFunction_x3<>::p2), 2.0);
  
  //2^3 = 8, dfdx: 3x^2, 3*2^2=12
  EQUAL(my_func(2), 8);
  EQUAL(interpolant->f_val(2), 8);
  EQUAL(interpolant->f_diff(2).f, 8);
  EQUAL(interpolant->f_diff(2).dfdx, 12);
  EQUAL(interpolant2->f_diff(2).f, 8);
  EQUAL(interpolant2->f_diff(2).dfdx, 12);
  
  interpolant->set_interval(-5,11);
  interpolant->set_size(8);
  double error;
  interpolant->interpolate_p1(error);
  
  EQUAL(interpolant->val(-7), -343);
  EQUAL(interpolant->val(-5), -125);
  EQUAL(interpolant->val(0), 0);
  EQUAL(interpolant->val(3), 27);
  EQUAL(interpolant->val(10), 1030);
  EQUAL(interpolant->val(11), 1331);
  
  EQUAL(interpolant->diff(-7).f, -343);
  EQUAL(interpolant->diff(-7).dfdx, 147);
  EQUAL(interpolant->diff(-4).dfdx, 51);
  EQUAL(interpolant->diff(3).dfdx, 27);
  EQUAL(interpolant->diff(5).dfdx, 75);
  EQUAL(interpolant->diff(9).dfdx, 243);
  EQUAL(interpolant->diff(10).dfdx, 303);
  EQUAL(interpolant->diff(11).dfdx, 363);
  
  delete interpolant;
  delete interpolant2;
}



