/*
 * mixed_test.cpp
 *
 *  Created on: Apr 8, 2019
 *      Author: jb
 */

#include <flow_gtest.hh>
#include <vector>
#include <iostream>
#include "tools/mixed.hh"

template <unsigned int dim, unsigned  int spacedim = 3>
class Mapping {
public:
    Mapping(int a, std::vector<int> b) : _a(a + dim + spacedim) {
    }

    Mapping(int a) : _a(a * dim * spacedim) {
    }

    int _a;
};


template <unsigned int dim>
class FE {
public:
    FE(int a, std::vector<int> b) : _a(a + dim + 10) {
    }

    FE(int a) : _a(a * dim * 10) {
    }

    int _a;
    int _b;
};


template <unsigned int dim>
class FE_XY : public FE<dim> {
public:
    FE_XY(int a, std::vector<int> b)
    : FE<dim>(a, b)
      { this->_b = 2*this->_a; }
};



TEST(Mixed, mixed) {
    // only dim parameter templates
    std::vector<int> vec = {1, 2};
    {
        // (int, vector<int>) constructor
        Mixed<FE> mixed_fe_a = Mixed<FE>( 3, vec);
        EXPECT_EQ(13, mixed_fe_a[Dim<0>{}]._a);
        EXPECT_EQ(14, mixed_fe_a[Dim<1>{}]._a);
        EXPECT_EQ(15, mixed_fe_a[2_d]._a);
        EXPECT_EQ(16, mixed_fe_a[3_d]._a);
    }
    {
        // (int) constructor
        Mixed<FE> mixed_fe_a = Mixed<FE>( 1 + 2);
        EXPECT_EQ(0, mixed_fe_a[Dim<0>{}]._a);
        EXPECT_EQ(30, mixed_fe_a[Dim<1>{}]._a);
        EXPECT_EQ(60, mixed_fe_a[2_d]._a);
        EXPECT_EQ(90, mixed_fe_a[3_d]._a);
    }


    {
        Mixed<FE_XY>  mixed_fe_xy( 3, vec);
        Mixed<FE> mixed_fe = mixed_fe_xy;
        EXPECT_EQ(13, mixed_fe[Dim<0>{}]._a);
        EXPECT_EQ(14, mixed_fe[Dim<1>{}]._a);
        EXPECT_EQ(15, mixed_fe[2_d]._a);
        EXPECT_EQ(16, mixed_fe[3_d]._a);
        EXPECT_EQ(26, mixed_fe[Dim<0>{}]._b);
        EXPECT_EQ(28, mixed_fe[Dim<1>{}]._b);
        EXPECT_EQ(30, mixed_fe[2_d]._b);
        EXPECT_EQ(32, mixed_fe[3_d]._b);
    }

// Compilation must failed on static assert (Non-convertible types!)
//    {
//        Mixed<FE>  mixed_fe(3);
//        Mixed<Mapping> mixed_map = mixed_fe;
//    }

    {
        // item constructor of Mixed
        FE<0> fe0(0);
        FE<1> fe1(1);
        FE<2> fe2(2);
        FE<3> fe3(3);
        Mixed<FE> mixed_fe_a = Mixed<FE>(fe0, fe1, fe2, fe3);
        EXPECT_EQ( 0, mixed_fe_a[Dim<0>{}]._a);
        EXPECT_EQ(10, mixed_fe_a[Dim<1>{}]._a);
        EXPECT_EQ(40, mixed_fe_a[2_d]._a);
        EXPECT_EQ(90, mixed_fe_a[3_d]._a);
    }

    // dim and spacedim templates
    {
        Mixed<Mapping> mixed_fe_a( 3, vec);
        EXPECT_EQ(6, mixed_fe_a[Dim<0>{}]._a);
        EXPECT_EQ(7, mixed_fe_a[Dim<1>{}]._a);
        EXPECT_EQ(8, mixed_fe_a[2_d]._a);
        EXPECT_EQ(9, mixed_fe_a[3_d]._a);
    }
    {
        Mixed<Mapping> mixed_fe_a( 3);
        EXPECT_EQ(0, mixed_fe_a[Dim<0>{}]._a);
        EXPECT_EQ(9, mixed_fe_a[Dim<1>{}]._a);
        EXPECT_EQ(18, mixed_fe_a[2_d]._a);
        EXPECT_EQ(27, mixed_fe_a[3_d]._a);
    }
//    {
//        MixedSpaceDim<Mapping> mixed_fe_a( 3, vec);
//        EXPECT_EQ(6, mixed_fe_a[Dim<0>{}]._a);
//        EXPECT_EQ(7, mixed_fe_a[Dim<1>{}]._a);
//        EXPECT_EQ(8, mixed_fe_a[2_d]._a);
//        EXPECT_EQ(9, mixed_fe_a[3_d]._a);
//    }
//    {
//        MixedSpaceDim<Mapping> mixed_fe_a( 3);
//        EXPECT_EQ(0, mixed_fe_a[Dim<0>{}]._a);
//        EXPECT_EQ(9, mixed_fe_a[Dim<1>{}]._a);
//        EXPECT_EQ(18, mixed_fe_a[2_d]._a);
//        EXPECT_EQ(27, mixed_fe_a[3_d]._a);
//    }

}

// Template specialization of Mixed<class T, 1> objects (doesn't contain dim=0)
TEST(Mixed, mixed_dim_1_3) {
    // only dim parameter templates
    std::vector<int> vec = {1, 2};
    {
        // (int, vector<int>) constructor
        Mixed<FE, 1> mixed_fe_a = Mixed<FE, 1>(3, vec);
        EXPECT_EQ(14, mixed_fe_a[Dim<1>{}]._a);
        EXPECT_EQ(15, mixed_fe_a[2_d]._a);
        EXPECT_EQ(16, mixed_fe_a[3_d]._a);
    }
    {
        // (int) constructor
        Mixed<FE, 1> mixed_fe_a = Mixed<FE, 1>( 1 + 2);
        EXPECT_EQ(30, mixed_fe_a[Dim<1>{}]._a);
        EXPECT_EQ(60, mixed_fe_a[2_d]._a);
        EXPECT_EQ(90, mixed_fe_a[3_d]._a);
    }


    {
        Mixed<FE_XY, 1>  mixed_fe_xy( 3, vec);
        Mixed<FE, 1> mixed_fe = mixed_fe_xy;
        EXPECT_EQ(14, mixed_fe[Dim<1>{}]._a);
        EXPECT_EQ(15, mixed_fe[2_d]._a);
        EXPECT_EQ(16, mixed_fe[3_d]._a);
        EXPECT_EQ(28, mixed_fe[Dim<1>{}]._b);
        EXPECT_EQ(30, mixed_fe[2_d]._b);
        EXPECT_EQ(32, mixed_fe[3_d]._b);
    }// */

// Compilation must failed on static assert (Non-convertible types!)
//    {
//        Mixed<FE, 1>  mixed_fe(3);
//        Mixed<Mapping, 1> mixed_map = mixed_fe;
//    }

    {
        // item constructor of Mixed
        FE<1> fe1(1);
        FE<2> fe2(2);
        FE<3> fe3(3);
        Mixed<FE, 1> mixed_fe_a = Mixed<FE, 1>(fe1, fe2, fe3);
        EXPECT_EQ(10, mixed_fe_a[Dim<1>{}]._a);
        EXPECT_EQ(40, mixed_fe_a[2_d]._a);
        EXPECT_EQ(90, mixed_fe_a[3_d]._a);
    }

    // dim and spacedim templates
    {
        Mixed<Mapping, 1> mixed_fe_a( 3, vec);
        EXPECT_EQ(7, mixed_fe_a[Dim<1>{}]._a);
        EXPECT_EQ(8, mixed_fe_a[2_d]._a);
        EXPECT_EQ(9, mixed_fe_a[3_d]._a);
    }
    {
        Mixed<Mapping, 1> mixed_fe_a( 3);
        EXPECT_EQ(9, mixed_fe_a[Dim<1>{}]._a);
        EXPECT_EQ(18, mixed_fe_a[2_d]._a);
        EXPECT_EQ(27, mixed_fe_a[3_d]._a);
    }

}

//void goo(int)
//{}
//
//void goo(int, float)
//{}
//
//template<typename... Args>
//void foo(Args&&... args)
//{std::cout << "arg_foo\n";
// goo(std::forward<Args>(args)...);
//}
//
//template < template<int dim> class TT>
//void foo( const MixedPtr<TT> &other)
//{std::cout << "dim_foo\n";}





TEST(MixedPtr, mixed_ptr) {
    // only dim parameter templates
    std::vector<int> vec = {1, 2};
    {
        MixedPtr<FE> mixed_fe = MixedPtr<FE>( 3, vec);
        EXPECT_EQ(13, mixed_fe[Dim<0>{}]->_a);
        EXPECT_EQ(14, mixed_fe[Dim<1>{}]->_a);
        EXPECT_EQ(15, mixed_fe[2_d]->_a);
        EXPECT_EQ(16, mixed_fe[3_d]->_a);
    }


    // dim and spacedim templates
    {
        MixedPtr<Mapping> mixed_fe( 3, vec);
        EXPECT_EQ(6, mixed_fe[Dim<0>{}]->_a);
        EXPECT_EQ(7, mixed_fe[Dim<1>{}]->_a);
        EXPECT_EQ(8, mixed_fe[2_d]->_a);
        EXPECT_EQ(9, mixed_fe[3_d]->_a);
    }
//    {
//        MixedSpaceDimPtr<Mapping> mixed_fe( 3, vec);
//        EXPECT_EQ(6, mixed_fe[Dim<0>{}]->_a);
//        EXPECT_EQ(7, mixed_fe[Dim<1>{}]->_a);
//        EXPECT_EQ(8, mixed_fe[2_d]->_a);
//        EXPECT_EQ(9, mixed_fe[3_d]->_a);
//    }

    // assign to base
    {
        MixedPtr<FE_XY>  fe_xy( 3, vec);
//        foo(fe_xy);
//        foo(3);
        auto mixed_fe = MixedPtr<FE>(fe_xy);
        EXPECT_EQ(13, mixed_fe[Dim<0>{}]->_a);
        EXPECT_EQ(14, mixed_fe[Dim<1>{}]->_a);
        EXPECT_EQ(15, mixed_fe[2_d]->_a);
        EXPECT_EQ(16, mixed_fe[3_d]->_a);
    }

    {
        // item constructor of MixedPtr
        std::shared_ptr< FE<0> > fe0 = std::make_shared< FE<0> >(0);
        std::shared_ptr< FE<1> > fe1 = std::make_shared< FE<1> >(1);
        std::shared_ptr< FE<2> > fe2 = std::make_shared< FE<2> >(2);
        std::shared_ptr< FE<3> > fe3 = std::make_shared< FE<3> >(3);
        MixedPtr<FE> mixed_fe = MixedPtr<FE>(fe0, fe1, fe2, fe3);
        EXPECT_EQ( 0, mixed_fe[Dim<0>{}]->_a);
        EXPECT_EQ(10, mixed_fe[Dim<1>{}]->_a);
        EXPECT_EQ(40, mixed_fe[2_d]->_a);
        EXPECT_EQ(90, mixed_fe[3_d]->_a);
    }

}


// Template specialization of MixedPtr<class T, 1> objects (doesn't contain dim=0)
TEST(MixedPtr, mixed_ptr_dim_1_3) {
    // only dim parameter templates
    std::vector<int> vec = {1, 2};
    {
        MixedPtr<FE, 1> mixed_fe = MixedPtr<FE, 1>( 3, vec);
        EXPECT_EQ(14, mixed_fe[Dim<1>{}]->_a);
        EXPECT_EQ(15, mixed_fe[2_d]->_a);
        EXPECT_EQ(16, mixed_fe[3_d]->_a);
    }


    // dim and spacedim templates
    {
        MixedPtr<Mapping, 1> mixed_fe( 3, vec);
        EXPECT_EQ(7, mixed_fe[Dim<1>{}]->_a);
        EXPECT_EQ(8, mixed_fe[2_d]->_a);
        EXPECT_EQ(9, mixed_fe[3_d]->_a);
    }
//    {
//        MixedSpaceDimPtr<Mapping> mixed_fe( 3, vec);
//        EXPECT_EQ(6, mixed_fe[Dim<0>{}]->_a);
//        EXPECT_EQ(7, mixed_fe[Dim<1>{}]->_a);
//        EXPECT_EQ(8, mixed_fe[2_d]->_a);
//        EXPECT_EQ(9, mixed_fe[3_d]->_a);
//    }

    // assign to base
    {
        MixedPtr<FE_XY, 1>  fe_xy( 3, vec);
//        foo(fe_xy);
//        foo(3);
        auto mixed_fe = MixedPtr<FE, 1>(fe_xy);
        EXPECT_EQ(14, mixed_fe[Dim<1>{}]->_a);
        EXPECT_EQ(15, mixed_fe[2_d]->_a);
        EXPECT_EQ(16, mixed_fe[3_d]->_a);
    }

    {
        // item constructor of MixedPtr
        std::shared_ptr< FE<1> > fe1 = std::make_shared< FE<1> >(1);
        std::shared_ptr< FE<2> > fe2 = std::make_shared< FE<2> >(2);
        std::shared_ptr< FE<3> > fe3 = std::make_shared< FE<3> >(3);
        MixedPtr<FE, 1> mixed_fe = MixedPtr<FE, 1>(fe1, fe2, fe3);
        EXPECT_EQ(10, mixed_fe[Dim<1>{}]->_a);
        EXPECT_EQ(40, mixed_fe[2_d]->_a);
        EXPECT_EQ(90, mixed_fe[3_d]->_a);
    }

}


//template<Dim dim>
//int foo(std::vector<float> &vec) {
//    return vec.size() * dim;
//}
//
//TEST(dim_switch, dim_switch) {
//    std::vector<int> vec = {1, 2};
//    for(int dim=0; dim<4; dim++) {
//        EXPECT_EQ(2+dim, dim_switch<foo>(dim, vec));
//    }
//}
