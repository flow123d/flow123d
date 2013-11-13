/*
 * time_governor_test.cpp
 *
 *  Created on: May 20, 2011
 *      Author: jb
 */

#include <gtest/gtest.h>

#include "interpolation/functors_impl.hh"
#include "interpolation/interpolant_impl.hh"

#define EQUAL(a,b) INPUT_CHECK( (a) == (b), #a": %f and "#b":%f differs\n",a,b);

//Example functor f(x)=x+2
template<class Type=double>
class Linear : public FunctorBase<Type>
{
public:
  ///Constructor.
  Linear(){}
  
  Type operator()(Type x)
  {
    return x+2;
  }    
};

//Example functor f(x)=x^2
template<class Type=double>
class Quadratic : public FunctorBase<Type>
{
public:
  ///Constructor.
  Quadratic(){}
  
  Type operator()(Type x)
  {
    return x*x;
  }    
};

//Example functor f(x)=x^3 + p1 
template<class Type=double>
class Cubic : public FunctorBase<Type>
{
public:
    
  typedef enum{ p1, p2, p3
  } Parameters;

  ///Constructor.
  Cubic(){}
  
  ///Constructor.
  Cubic(double pp1, double pp2, double pp3)
  {
    this->set_param(p1, pp1);
    this->set_param(p2, pp2);
    this->set_param(p3, pp3);
  }
  
  Type operator()(Type x)
  {
    return x*x*x + this->param(p1);
  }    
};


//Example functor f(x)=x^2 + y^2 - 4 
template<class Type=double>
class Circle : public IFunctorBase<Type>
{
public:
    
  typedef enum{ radius, c_x, c_y
  } Parameters;

  ///Constructor.
  Circle(double prad, double pc_x, double pc_y)
  {
    this->set_param(radius, prad);
    this->set_param(c_x, pc_x);
    this->set_param(c_y, pc_y);
  }
  
  Type operator()(Type x,Type y)
  {
    return std::pow(x-this->param(c_x),2) + std::pow(y-this->param(c_y),2) - std::pow(this->param(radius),2);
  }    
};


TEST (Functors, functors)
{
  //Explicit functor and setting parameters
  Cubic<double> my_func;                //x^3

  //setting parameters
  my_func.set_param(Cubic<double>::p1,1.0);
  my_func.set_param(Cubic<>::p2,2.0);
  my_func.set_param(Cubic<>::p3,3.0);
  EQUAL(my_func.param(2), 3.0);
  EQUAL(my_func.param(Cubic<>::p2), 2.0);
  EQUAL(my_func(2.0), 9.0);
  
  //Explicit functor and setting parameters in constructor
  Cubic<> my_func2(3.0,4.0,5.0);        //x^3+2 (parameter p1=2.0)
  EQUAL(my_func2.param(2), 5.0);
  EQUAL(my_func2.param(Cubic<>::p1), 3.0);
  EQUAL(my_func2(2.0), 11.0);
  
  
  //copying functor and parameters
  Cubic<double> my_func3;
  my_func3.set_param_from_func(&my_func2);
  
  EQUAL(my_func3.param(2), 5.0);
  EQUAL(my_func3.param(Cubic<>::p1), 3.0);
  
  
  //Implicit functor
  Circle<double> circle_func(4.0, 1.0, 2.0);    //(x-1)^2 + (y-2)^2 - 4^2 = 0
  
  EQUAL(circle_func.param(Circle<>::radius), 4.0);
  EQUAL(circle_func.param(Circle<>::c_x), 1.0);
  EQUAL(circle_func(5.0,2.0), 0.0);
  
  /*
  InterpolantImplicit* interpolant = new InterpolantImplicit(&circle_func);
  
  interpolant->fix_variable(InterpolantImplicit::fix_x,2.0);
  
  EQUAL(interpolant->f_val(3.0),-3.0);           //2^2+3^2-4^2 = 4+9-16=-3
  //*/
}


/**
 * Test for interpolant of x^3
 */
TEST (Functors, make_interpolation)
{
  Cubic<double> my_func;                //x^3
  Cubic<> my_func2;                     //x^3+p1 (parameter p1=2.0)
  
  my_func.set_param(Cubic<double>::p1,0.0);
  my_func2.set_param(Cubic<double>::p1,2.0);
  
  Interpolant* interpolant = new Interpolant(&my_func);
  Interpolant* interpolant2 = new Interpolant();
  interpolant2->set_functor(&my_func2);
  //extrapolation set by default to functor

  //statistics - zero at start
  EQUAL(interpolant->statistics().total_calls, 0);
  EQUAL(interpolant->statistics().interval_miss_a, 0);
  EQUAL(interpolant->statistics().interval_miss_b, 0);
  EQUAL(interpolant->statistics().min, 0);
  EQUAL(interpolant->statistics().max, 0);
  
  //2^3 = 8, dfdx: 3x^2, 3*2^2=12
  EQUAL(my_func(2), 8);
  EQUAL(interpolant->f_val(2), 8);
  EQUAL(interpolant->f_diff(2).first, 8);
  EQUAL(interpolant->f_diff(2).second, 12);
  
  //2nd 3rd derivate: 6*x, 6
  EQUAL(interpolant->f_diffn(4,2), 24);
  EQUAL(interpolant->f_diffn(4,3), 6);
  
  EQUAL(interpolant2->f_diff(2).first, 10);
  EQUAL(interpolant2->f_diff(2).second, 12);
  
  
  interpolant->set_interval(-5,11);
  interpolant->set_size(8);
  interpolant->interpolate();
  
  //DBGMSG("Error of interpolation: %f\n", interpolant->error());
  
  EQUAL(interpolant->val(-7), -343);    //out of interval
  EQUAL(interpolant->statistics().min, -7.0);
  EQUAL(interpolant->statistics().max, 11);
  
  EQUAL(interpolant->val(-5), -125);
  EQUAL(interpolant->val(0), 0);
  EQUAL(interpolant->val(3), 27);
  EQUAL(interpolant->val(10), 1030);
  EQUAL(interpolant->val(11), 1331);
  EQUAL(interpolant->val(15), 3375);    //out of interval
  
  EQUAL(interpolant->diff(-7).first, -343); //out of interval
  EQUAL(interpolant->diff(-7).second, 147);       //out of interval
  EQUAL(interpolant->diff(-4).second, 51);
  EQUAL(interpolant->diff(3).second, 27);
  EQUAL(interpolant->diff(5).second, 75);
  EQUAL(interpolant->diff(9).second, 243);
  EQUAL(interpolant->diff(10).second, 303);
  EQUAL(interpolant->diff(11).second, 363);
  
  //statistics
  EQUAL(interpolant->statistics().total_calls, 15);
  EQUAL(interpolant->statistics().interval_miss_a, 3);
  EQUAL(interpolant->statistics().interval_miss_b, 1);
  EQUAL(interpolant->statistics().min, -7.0);
  EQUAL(interpolant->statistics().max, 15.0);
  
  //extrapolation
  //functor type has been tested as default before
  interpolant->set_extrapolation(Extrapolation::constant);
  EQUAL(interpolant->val(-10), -125);
  EQUAL(interpolant->val(15), 1331);
  EQUAL(interpolant->diff(-10).first, -125);
  EQUAL(interpolant->diff(15).first, 1331);
  EQUAL(interpolant->diff(-10).second, 75);
  EQUAL(interpolant->diff(15).second, 363);
  interpolant->set_extrapolation(Extrapolation::linear);
  EQUAL(interpolant->val(-10), -370);           //-125 + 49*(-10-(-5))
  EQUAL(interpolant->val(15), 2535);            //729 + 301*(15-9)
  EQUAL(interpolant->diff(-10).first, -370);    
  EQUAL(interpolant->diff(15).first, 2535);
  EQUAL(interpolant->diff(-10).second, 195);    //75 + (-24)*(-10-(-5))
  EQUAL(interpolant->diff(15).second, 603);     //243 + 60*(15-9)
  
  //increasing missed interval evaluation
  //this should remake the table
  for(unsigned int i=0; i < (STATISTIC_CHECK+1); i++)
  {
    //DBGMSG("i: %d\n",i);
    interpolant->val(-20);
    interpolant->val(15);
  }
  
  interpolant->check_stats_and_reinterpolate();
  //interpolant->check_stats_and_reinterpolate(0.5);
  
  //statistics
  //DBGMSG("bound_a: %f\n",interpolant->bound_a());
  //DBGMSG("bound_b: %f\n",interpolant->bound_b());
  //DBGMSG("total: %d\n",interpolant->statistics().total_calls);
  EQUAL(interpolant->bound_a(),-20);
  EQUAL(interpolant->bound_b(),15);
  
  EQUAL(interpolant->statistics().total_calls, 0);
  EQUAL(interpolant->statistics().interval_miss_a, 0);
  EQUAL(interpolant->statistics().interval_miss_b, 0);
  EQUAL(interpolant->statistics().min, -20.0);
  EQUAL(interpolant->statistics().max, 15.0);

  delete interpolant;
  delete interpolant2;
}




TEST (Functors, interpolation_error)
{
  //note: when computing integral of the difference don't forget
  //that the variable is moved: x'=(x-x[i])
  
  Linear<double> lin_func;                //x+2
  Quadratic<double> quad_func;            //x^2
  Cubic<double> cubic_func;               //x^3
  
  Interpolant* interpolant = new Interpolant(&lin_func);
  interpolant->set_interval(0,10);
  interpolant->set_size(2);
  
  EQUAL(interpolant->error(), -1);      //error not computed yet
  
  interpolant->set_functor<Linear, double>(&lin_func);
  interpolant->interpolate();
  //linear function is interpolated by linear aproximation accurately
  //DBGMSG("Error of interpolation: %.64f\n", interpolant->error());
  EQUAL(interpolant->error(), 0);
  
    /*
  interpolant->set_functor<Quadratic, double>(&quad_func);
  interpolant->interpolate();
  //http://www.numberempire.com/definiteintegralcalculator.php
  // \int_0^5 (x^2-5*x)^2+(2*x-2*x)^2 = 104.1666666666666
  // sqrt = 10,206207262
  // \int_5^10 (x^2-15*(x-5)-25)^2+(2*x-2*(x-5)-10)^2 = 104.1666666666666
  // sqrt = 10,206207262
  //DBGMSG("Error of interpolation: %.64f\n", interpolant->error());
  EQUAL(interpolant->error(), 2.041241452319315197172500120359472930431365966796875);
  
  
  cubic_func.set_param(Cubic<double>::p1,0.0);
  interpolant->set_functor<Cubic, double>(&cubic_func);
  interpolant->interpolate();
  //http://www.numberempire.com/definiteintegralcalculator.php
  // \int_0^5 (x^3-25*x)^2+(3*x^2-15*x)^2 = 6889.880952380952
  // sqrt = 83,005306772
  // \int_5^10 (x^3-175*(x-5)-125)^2+(3*x^2-45*(x-5)-75)^2 = 53764.88095238095
  // sqrt = 231,872553254
  //DBGMSG("Error of interpolation: %.64f\n", interpolant->error());
  EQUAL(interpolant->error(), 31.48778600260849458436496206559240818023681640625);
  
  interpolant->set_size_automatic(5e-3,5);
  interpolant->interpolate();
  //DBGMSG("Error of interpolation: %.64f\n", interpolant->error());
  */
  delete interpolant;
}


