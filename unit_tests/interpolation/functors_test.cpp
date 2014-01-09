/*
 * functors_test.cc
 *
 *  Created on: November, 2013
 *      Author: pe
 */

#include <ctime>
#include <cstdlib>
#include <sstream>

#include <flow_gtest.hh>

#include "interpolation/functors_impl.hh"
#include "interpolation/interpolant_impl.hh"

#define EQUAL(a,b) INPUT_CHECK( (a) == (b), #a": %f and "#b":%f differs\n",a,b);

//*********************************** FUNCTOR EXAMPLES *********************************

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


//Example functor f(x)=x^3
template<class Type=double>
class Power3 : public FunctorBase<Type>
{
public:
    
  typedef enum{ p1, p2, p3
  } Parameters;

  ///Constructor.
  Power3(){}
  
  Type operator()(Type x)
  {
    return pow(x,3);
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
  Circle(){}
  
  ///Constructor.
  Circle(double prad, double pc_x, double pc_y)
  {
    this->set_param(radius, prad);
    this->set_param(c_x, pc_x);
    this->set_param(c_y, pc_y);
  }
  
  Type operator()(Type x,Type y)
  {
    //return std::pow(x-this->param(c_x),2) + std::pow(y-this->param(c_y),2) - std::pow(this->param(radius),2);
    return pow(x-this->param(c_x),2) + pow(y-this->param(c_y),2) - pow(this->param(radius),2);
  }    
};


//*********************************** UNIT TESTS *********************************

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
}



TEST (Functors, make_interpolation)
{
  bool interpolate_derivative = true;
  Cubic<double> my_func;                //x^3
  Cubic<> my_func2;                     //x^3+p1 (parameter p1=2.0)
  
  //test interpolant2-----------
  my_func2.set_param(Cubic<double>::p1,2.0);
  Interpolant interpolant2;
  interpolant2.set_functor(&my_func2, false);
  interpolant2.set_interval(0, 10);
  //interpolant2.set_size(10);
  interpolant2.set_size_automatic(0.01, 10, 1e5);
  EQUAL(interpolant2.f_diff(2).first, 10);
  EQUAL(interpolant2.f_diff(2).second, 12);
  interpolant2.interpolate();
  EXPECT_TRUE(interpolant2.error() < 3.891645e-3);
  EXPECT_TRUE(interpolant2.error() > 3.891643e-3);
  //-----------------------
  
  
  my_func.set_param(Cubic<double>::p1,0.0);
  Interpolant* interpolant = new Interpolant(&my_func, interpolate_derivative);
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
  EXPECT_EQ(interpolant->f_diffn(4,2), 24);
  EXPECT_EQ(interpolant->f_diffn(4,3), 6);
  
  interpolant->set_interval(-5,11);
  interpolant->set_size(8);
  //interpolant->set_size_automatic(0.1,10,1e5);
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
  for(unsigned int i=0; i < 100; i++)
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
}
//*/

TEST(Functors, make_implicit)
{
  Circle<double> circle_func(4.0, 1.0, 2.0);    //(x-1)^2 + (y-2)^2 - 4^2 = 0
  
  
  InterpolantImplicit interpolant(&circle_func, true);
  
  interpolant.fix_variable(IFixVariable::fix_x,2.0);
  
  EQUAL(interpolant.f_val(3.0),-14.0);           //(2-1)^2+(3-2)^2-4^2 = 1+1-16=-14
  
  interpolant.set_interval(-4,4);
  interpolant.set_size(10);
  //interpolant.interpolate();
}


TEST (Functors, interpolation_error)
{
  //note: when computing integral of the difference don't forget
  //that the variable is moved: x'=(x-x[i])
  
  Linear<double> lin_func;                //x+2
  Quadratic<double> quad_func;            //x^2
  Cubic<double> cubic_func;               //x^3
  
  Interpolant* interpolant = new Interpolant(&lin_func);
  interpolant->set_interval(1,10);
  interpolant->set_size(10);
  
  EQUAL(interpolant->error(), -1);      //error not computed yet
  
  //interpolant->set_functor<Linear, double>(&lin_func);
  interpolant->interpolate();
  //linear function is interpolated by linear aproximation accurately
  //DBGMSG("Error of interpolation: %.64f\n", interpolant->error());
  //cout << "Error of interpolation: " << interpolant->error() << endl;
  //EQUAL(interpolant->error(), 0);
  EXPECT_TRUE(interpolant->error() < 1e-15);
  
  interpolant->set_functor<Quadratic, double>(&quad_func);
  interpolant->interpolate();
  EXPECT_TRUE(interpolant->error() < 9.631392E-02);
  EXPECT_TRUE(interpolant->error() > 9.631390E-02);
  
  
  
  cubic_func.set_param(Cubic<double>::p1,0.0);
  interpolant->set_functor<Cubic, double>(&cubic_func);
  interpolant->set_interval(1,10);
  interpolant->set_norm(ErrorNorm::l2);
  interpolant->interpolate();
  
  
  delete interpolant;
}


//*********************************** FIXTURE TESTS - BENCHMARKS *********************************

class InterpolantCreateTest : public ::testing::Test{
protected:
  virtual void SetUp()
  {
    size = 3*1000*1000;
    a = 1;
    b = 10;

    interpolant.set_functor(&cubic_func, false);
    interpolant.set_interval(a,b);
  }
  
  Power3<double> cubic_func;
  //virtual void TearDown() {}
  double a, b, step;
  unsigned int size;
  Interpolant interpolant;
};

TEST_F(InterpolantCreateTest, Fix)
{
  double step = (b-a)/size;
  
  vector<double> x_vec(size+1);
  vector<double> f_vec(size+1);
  vector<double> p_vec(size);
  
  double temp_x = a;
  //filling the vector x and f
  for(unsigned int i = 0; i != size; i++)
  {
    temp_x = a + step*i;
    x_vec[i] = temp_x;
    f_vec[i] = cubic_func(temp_x);
  }
  //finish the interval
  x_vec[size] = b;
  f_vec[size] = cubic_func(b);

  for(unsigned int i = 0; i != size; i++)
  {
    p_vec[i] = (f_vec[i+1] - f_vec[i]) / (x_vec[i+1] - x_vec[i]);
  }

  //to this point 254ms
   
  //error computation in the middle nodes
  //we do not create vector of middle nodes here
  //so that can make some difference
  double temp_a = a+step/2;
  double tot_err, value, int_value;
  for(unsigned int i = 0; i != size; i++)
  {
    temp_x = temp_a + step*i;
    value = cubic_func(temp_x);
    unsigned int ii = floor((temp_x - a) / step);
    int_value = (temp_x - x_vec[ii])*p_vec[ii] + f_vec[ii];
    tot_err = std::max( tot_err,
                        std::abs(value - int_value) / (std::abs(value)+1e-10)
                      );
  }
  //to this point 467ms
};


TEST_F(InterpolantCreateTest, Interpolate)
{
  interpolant.set_size(size);
  interpolant.interpolate();
};

TEST_F(InterpolantCreateTest, Automatic)
{
  interpolant.set_size_automatic(1e-11, 11, size);
  interpolant.interpolate();
};







class InterpolantTest : public ::testing::Test{
protected:
  virtual void SetUp()
  {
    unsigned int size = 10000;
    a = 1;
    b = 10;
    
    n = 10*1000*1000;
    
    step = (b-a)/n;
    
    interpolant.set_functor(&cubic_func, false);
    interpolant.set_interval(a,b);
    interpolant.set_size(size);
    interpolant.interpolate();
  }
  
  Power3<double> cubic_func;
  //virtual void TearDown() {}
  double a, b, step;
  unsigned int n;
  Interpolant interpolant;
};


// 3.12.2013, Pavel Exner
// 
// Test results:
// [       OK ] InterpolantTest.FunctorEval (679 ms)
// [       OK ] InterpolantTest.InterpolantEval_val (65 ms)
// [       OK ] InterpolantTest.InterpolantEval_val_p1 (196 ms)
// [       OK ] InterpolantTest.InterpolantEval_val_test (194 ms)
// 
// For unknown reason the val() function is faster than val_p1() which is called inside val().
// It also behave differently when only some of these 4 tests are running.
// 
// Test results of valgrind run:
// > valgrind --tool=callgrind --dump-instr=yes ./functors_test_bin
// 
// [       OK ] InterpolantTest.FunctorEval (19593 ms)
// [       OK ] InterpolantTest.InterpolantEval_val (3640 ms)
// [       OK ] InterpolantTest.InterpolantEval_val_p1 (2470 ms)
// [       OK ] InterpolantTest.InterpolantEval_val_test (2496 ms)
// 
// It is also strange that the val() function is slower than the others when using Valgrind.
// 
// Conclusion:
// On the other hand, the interpolation is still at least 3 times faster than 
// the functor evaluation (std::pow) and that can be satisfying for now...


// 5.12.2013, Pavel Exner
//
// Test results:
// [       OK ] InterpolantTest.FunctorEval (659 ms)
// [       OK ] InterpolantTest.InterpolantEval_val (165 ms)
// [       OK ] InterpolantTest.InterpolantEval_val_p1 (164 ms)
// [       OK ] InterpolantTest.InterpolantEval_val_test (162 ms)
//
// After removing unnecessary node and coeficient arrays:
// These results are satisfying:
//  + the function val_p1 itself is by 14% faster than before
//  + we use only one array to keep the function values
//  + all three functions take the same time which means
//    that the collecting statistics and going through several IF commands
//    does not affect the evaluation time


TEST_F(InterpolantTest, FunctorEval){
  double res = 0;
  for(unsigned int i= 0; i < n; i++)
    {
      res += cubic_func(a+step*i);
    }
  cout << "res = " << res << endl;
}


TEST_F(InterpolantTest, InterpolantEval_val){
  double res = 0;
  for(unsigned int i= 0; i < n; i++)
    {
      res += interpolant.val(a+step*i);
    }
  cout << "res = " << res << endl;
}


TEST_F(InterpolantTest, InterpolantEval_val_p1){
  double res = 0;
  for(unsigned int i= 0; i < n; i++)
    {
      res += interpolant.val_p1(a+step*i);
    }
  cout << "res = " << res << endl;
}


TEST_F(InterpolantTest, InterpolantEval_val_test){
  double res = 0;
  for(unsigned int i= 0; i < n; i++)
    {
      res += interpolant.val_test(a+step*i);
    }
  cout << "res = " << res << endl;
}



//***************************   HYDROLOGY FUNCTIONS   *****************************

template<class Type=double>
class FK : public FunctorBase<Type>   //FK - hydraulic conductivity function
{
  private:
    double Bpar,PPar,n,Qr,Qs,
           Qa,Qm,Qk,Alfa,Ks,Kk,
           m,Hr,Hk,Hs,
           C1Qee,C2Qee,Qeer,Qees,Qeek;
        
  public:
    FK()
    {
      Bpar = 0.5,
      PPar = 2;
      
      n = 1.111,
      Qr = 0.001,
      Qs = 0.436,
      Qa = 0.001,
      Qm = 0.439,
      Qk = 0.436,
      Alfa = 0.733,
      Ks = 0.0505,
      Kk = 0.0505,
        
                
      m = 1-1/n;            
      C1Qee = 1/(Qm - Qa);
      C2Qee = -Qa/(Qm - Qa);
      Qeer = max(C1Qee*Qr + C2Qee,1E-3);
      Qees = min(C1Qee*Qs + C2Qee,1 - 1E-15);
      Qeek = min(C1Qee*Qk + C2Qee,Qees);
      Hr = -1/Alfa*pow(pow(Qeer,-1/m)-1,1/n);
      Hs = -1/Alfa*pow(pow(Qees,-1/m)-1,1/n); 
      Hk = -1/Alfa*pow(pow(Qeek,-1/m)-1,1/n); 
    }
    
    double get_Hr()
    { return Hr; }
    
    double get_Hk()
    { return Hk; }
    
    double get_Hs()
    { return Hs; }
    
    
    Type operator()(Type h)
    {   
      Type Kr,FFQr,FFQ,FFQk,Qee,Qe,Qek,C1Qe,C2Qe,Q;
      
      if (h <= Hr) return Ks*(1E-9);
      else if(h < Hk)
      {
            Q = Qa + (Qm - Qa)*pow((1 + pow(-Alfa*h,n)),-m);
            Qee = C1Qee*Q + C2Qee;
            FFQr = pow(1 - pow(Qeer,1/m),m);
            FFQ = pow(1 - pow(Qee,1/m),m);
            FFQk = pow(1 - pow(Qeek,1/m),m);
            C1Qe = 1/(Qs - Qr);
            C2Qe = -Qr/(Qs - Qr);
            Qe = C1Qe*Q + C2Qe;
            Qek = C1Qe*Qk + C2Qe;
            Kr = pow(Qe/Qek,Bpar)*pow((FFQr - FFQ)/(FFQr - FFQk),PPar) * Kk/Ks;
            return max<Type>(Ks*Kr,Ks*(1E-9));  
      }
      else if(h <= Hs)
      {
           Kr = (1-Kk/Ks)/(Hs-Hk)*(h-Hs) + 1;
           return Ks*Kr;
      }        
      else return Ks;
   }               
        
};


template<class Type=double>
class FQ : public FunctorBase<Type>            //FQ - water saturation function
{
  private:
    double n,Qr,Qs,Qa,Qm,Alfa,
           m,C1Qee,C2Qee,Qeer,Qees,Hr,Hs;
    
  public:
    FQ()
    {
        n = 1.111,
        Qr = 0.001,
        Qs = 0.436,
        Qa = 0.001,
        Qm = 0.439,
        Alfa = 0.733,
    
        m = 1 - 1/n;
        C1Qee = 1/(Qm - Qa);
        C2Qee = -Qa/(Qm - Qa);
        Qeer = max(C1Qee*Qr + C2Qee,1E-3);
        Qees = min(C1Qee*Qs + C2Qee,1-1E-15);
        Hr = -1/Alfa*pow(pow(Qeer,-1/m)-1,1/n);
        Hs = -1/Alfa*pow(pow(Qees,-1/m)-1,1/n); 
    }
    
    double get_Hr()
    { return Hr; }
    
    double get_Hs()
    { return Hs; }
        
    Type operator()( Type h )
    {
      Type Qee;
      if(h <= Hr) return Qr;
      else if(h < Hs)
      {
        Qee = pow(1+pow(-Alfa*h,n),-m);
        return Qa + (Qm - Qa)*Qee;
      }
      else return Qs;
    }
};


class HydrologyTest : public ::testing::Test{
protected:
  virtual void SetUp()
  {
    unsigned int size = 10000;
    a = -10.0;
    b = fk.get_Hs()-1e-10;        //-0.127;
    
    n = 6*1000*1000;
    
    step = (b-a)/n;
    
    interpolant_fk.set_functor(&fk);
    interpolant_fq.set_functor(&fq);
    interpolant_fk.set_interval(a,b);
    interpolant_fq.set_interval(a,b);
    interpolant_fk.set_size(size);
    interpolant_fq.set_size(size);
  }
  
  FK<double> fk;
  FQ<double> fq;
  //virtual void TearDown() {}
  double a, b, step;
  unsigned int n;
  Interpolant interpolant_fk;
  Interpolant interpolant_fq;
};


// 5.12.2013, Pavel Exner
// 
// Results of following tests with n=6*1000*1000 and size=10000.
// [       OK ] HydrologyTest.FKEval (3603 ms)
// [       OK ] HydrologyTest.FQEval (824 ms)
// [       OK ] HydrologyTest.Interpolate_FK (113 ms)
// [       OK ] HydrologyTest.Interpolate_FQ (103 ms)


TEST_F(HydrologyTest, FKEval){
  double res = 0;
  for(unsigned int i= 0; i < n; i++)
    {
      res += fk(a+step*i);
    }
  cout << "res = " << res << endl;
}

TEST_F(HydrologyTest, FQEval){
  double res = 0;
  for(unsigned int i= 0; i < n; i++)
    {
      res += fq(a+step*i);
    }
  cout << "res = " << res << endl;
}

TEST_F(HydrologyTest, Interpolate_FK){
  interpolant_fk.interpolate();
  
  double res = 0;
  for(unsigned int i= 0; i < n; i++)
    {
      res += interpolant_fk.val(a+step*i);
    }
  cout << "res = " << res << endl;
}

TEST_F(HydrologyTest, Interpolate_FQ){
  interpolant_fq.interpolate();
  
  double res = 0;
  for(unsigned int i= 0; i < n; i++)
    {
      res += interpolant_fq.val(a+step*i);
    }
  cout << "res = " << res << endl;
}