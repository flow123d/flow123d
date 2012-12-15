/*
 * field_value_test.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#include <gtest/gtest.h>
#include <fields/field_values.hh>

#include <iostream>
using namespace std;

TEST(FieldValue_, all) {

    {
        typedef FieldValue_<2,2, double> T;
        T::return_type x_val;
        T val(x_val);
        val(0,0)=3.0; val(0,1)=4.0;
        EXPECT_EQ(3.0,  val(0,0) );
        EXPECT_EQ(4.0,  val(0,1) );

        T::return_type r_val = val;
        EXPECT_EQ(3.0,  r_val(0,0) );
        EXPECT_EQ(4.0,  r_val(0,1) );

    }

}


/**
 * Speed results:
 * debug (-g -O0 -NODEBUG) (100 M steps):
 * interface: 1747ms
 * direct   :  361ms
 *
 * optimized -O3 (100 M steps):
 * interface: 123ms
 * direct   : 121ms
 */

TEST(FieldValue_, speed_test_interface) {

   typedef FieldValue_<1,1, double> T;
   double r_val;


   for(int step=0;step < 100000000; step++) {
       T val(r_val);

       for(int row=0;row< val.n_cols(); ++row)
           for(int col=0;col< val.n_rows(); ++col)
               val(row,col)+=step;
   }
   cout << r_val << endl;

}

TEST(FieldValue_, speed_test_direct) {

   double val;

   for(int step=0;step < 100000000; step++) {
       val+=step;
   }
   cout << val << endl;
}
