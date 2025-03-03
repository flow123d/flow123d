/*
 * backtrace_test.cpp
 *
 *  Created on: Jan 28, 2013
 *      Author: jb
 */



#include <flow_gtest.hh>

// in the third_party/FADBAD++ dir, namespace "fadbad"
#include "fadbad.h"
#include "badiff.h"
#include "fadiff.h"

#include "flow/soil_models.hh"

using namespace std;

template <class Model>
void water_content(Model m, double head, double wc, double cap) {
    fadbad::B<double> x_phead(head);
    fadbad::B<double> evaluated = m.water_content_diff(x_phead);
    evaluated.diff(0,1);
    double water_content = evaluated.x();
    double capacity = x_phead.d(0);
    cout << "p: " << head << " wc: " << water_content << " cap: " << capacity << endl << endl;
    EXPECT_NEAR(wc, water_content, 1e-5);
    EXPECT_NEAR(cap, capacity, 1e-5);

}

template <class Model>
void conductivity(Model m, double head, double cond, double d_cond) {
    //double conductivity_ = m.conductivity(head);
    //double d_conductivity = 0;

    fadbad::B<double> x_phead(head);
    fadbad::B<double> evaluated( m.conductivity_diff(x_phead) );
    evaluated.diff(0,1);
    double conductivity_ = evaluated.val();
    double d_conductivity = x_phead.d(0);
    cout << "p: " << head << " c: " << conductivity_ << " dc: " << d_conductivity << endl << endl;
    EXPECT_NEAR(cond, conductivity_, 1e-14);
    EXPECT_NEAR(d_cond, d_conductivity, 1e-14);
}

template <class Model>
void check(Model m, double head, double wc, double cap, double cond, double d_cond) {
    water_content(m, head, wc, cap);
    conductivity(m, head, cond, d_cond);
}


TEST(soil_model_VanGenuchten, all) {
    SoilModel_VanGenuchten soil_model;
    // bentonit
    SoilData soil_data;
    soil_data.n = 1.24;
    soil_data.alpha = 0.005;
    soil_data.Qr = 0.04;
    soil_data.Qs = 0.42;
    soil_data.Ks = 1.8e-10;
    soil_data.cut_fraction = 0.999;

    soil_model.reset(soil_data);


      check(soil_model, 1.00000000e+00, 4.20000000e-01, 0.00000000e+00, 1.80000000e-10, 0.00000000e+00);
       check(soil_model, 0.00000000e+00, 4.20000000e-01, 0.00000000e+00, 1.80000000e-10, 0.00000000e+00);
       check(soil_model, -1.00000000e+00, 4.20000000e-01, 0.00000000e+00, 1.80000000e-10, 0.00000000e+00);
       check(soil_model, -3.00000000e+00, 4.19978640e-01, 1.65512662e-04, 1.77846320e-10, 1.62768465e-11);
       check(soil_model, -1.00000000e+01, 4.18612382e-01, 2.16110151e-04, 1.16624573e-10, 5.18280353e-12);
       check(soil_model, -1.00000000e+02, 3.95257731e-01, 2.53605674e-04, 1.86364343e-11, 2.44082904e-13);
       check(soil_model, -1.00000000e+03, 2.92204588e-01, 5.32865027e-05, 2.13029689e-13, 5.08405076e-16);
       check(soil_model, -1.00000000e+04, 1.88528718e-01, 3.53702521e-06, 6.25249530e-16, 9.30855868e-19);


}



TEST(soil_model_Irmay, all) {
    SoilModel_Irmay soil_model;
    // bentonit
    SoilData soil_data;
    soil_data.n = 1.24;
    soil_data.alpha = 0.005;
    soil_data.Qr = 0.04;
    soil_data.Qs = 0.42;
    soil_data.Ks = 1.8e-10;
    soil_data.cut_fraction = 0.999;

    soil_model.reset(soil_data);


     check(soil_model, 1.00000000e+00, 4.20000000e-01, 0.00000000e+00, 1.80000000e-10, 0.00000000e+00);
      check(soil_model, 0.00000000e+00, 4.20000000e-01, 0.00000000e+00, 1.80000000e-10, 0.00000000e+00);
      check(soil_model, -1.00000000e+00, 4.20000000e-01, 0.00000000e+00, 1.80000000e-10, 0.00000000e+00);
      check(soil_model, -3.00000000e+00, 4.19978640e-01, 1.65512662e-04, 1.79430279e-10, 2.34454005e-13);
      check(soil_model, -1.00000000e+01, 4.18612382e-01, 2.16110151e-04, 1.77501741e-10, 3.04247369e-13);
      check(soil_model, -1.00000000e+02, 3.95257731e-01, 2.53605674e-04, 1.46638764e-10, 3.13553151e-13);
      check(soil_model, -1.00000000e+03, 2.92204588e-01, 5.32865027e-05, 5.24659003e-11, 3.34749654e-14);
      check(soil_model, -1.00000000e+04, 1.88528718e-01, 3.53702521e-06, 1.07164139e-11, 7.75481824e-16);
}
