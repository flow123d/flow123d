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
};


template <unsigned int dim>
class FE_XY : public FE<dim> {
public:
    FE_XY(int a, std::vector<int> b)
    : FE<dim>(a, b)
      {}
};



TEST(Mixed, mixed) {
    // only dim parameter templates
    std::vector<int> vec = {1, 2};
    {
        // (int, vector<int>) constructor
        Mixed<FE> mixed_fe_a = Mixed<FE>( 3, vec);
        EXPECT_EQ(13, mixed_fe_a.get<0>()._a);
        EXPECT_EQ(14, mixed_fe_a.get<1>()._a);
        EXPECT_EQ(15, mixed_fe_a.get<2>()._a);
        EXPECT_EQ(16, mixed_fe_a.get<3>()._a);
    }
    {
        // (int) constructor
        Mixed<FE> mixed_fe_a = Mixed<FE>( 1 + 2);
        EXPECT_EQ(0, mixed_fe_a.get<0>()._a);
        EXPECT_EQ(30, mixed_fe_a.get<1>()._a);
        EXPECT_EQ(60, mixed_fe_a.get<2>()._a);
        EXPECT_EQ(90, mixed_fe_a.get<3>()._a);
    }


    {
        Mixed<FE_XY>  mixed_fe_xy( 3, vec);
        Mixed<FE> mixed_fe = Mixed<FE>::cast_to_parent_template<FE_XY>(mixed_fe_xy);
        EXPECT_EQ(13, mixed_fe.get<0>()._a);
        EXPECT_EQ(14, mixed_fe.get<1>()._a);
        EXPECT_EQ(15, mixed_fe.get<2>()._a);
        EXPECT_EQ(16, mixed_fe.get<3>()._a);
    }

    /*{
        Mixed<FE_XY>  mixed_fe_xy( 3, vec);
        Mixed<FE> mixed_fe = mixed_fe_xy.template operator()<FE>();
        EXPECT_EQ(13, mixed_fe.get<0>()._a);
        EXPECT_EQ(14, mixed_fe.get<1>()._a);
        EXPECT_EQ(15, mixed_fe.get<2>()._a);
        EXPECT_EQ(16, mixed_fe.get<3>()._a);
    } // */

    // dim and spacedim templates
    {
        Mixed<Mapping> mixed_fe_a( 3, vec);
        EXPECT_EQ(6, mixed_fe_a.get<0>()._a);
        EXPECT_EQ(7, mixed_fe_a.get<1>()._a);
        EXPECT_EQ(8, mixed_fe_a.get<2>()._a);
        EXPECT_EQ(9, mixed_fe_a.get<3>()._a);
    }
    {
        Mixed<Mapping> mixed_fe_a( 3);
        EXPECT_EQ(0, mixed_fe_a.get<0>()._a);
        EXPECT_EQ(9, mixed_fe_a.get<1>()._a);
        EXPECT_EQ(18, mixed_fe_a.get<2>()._a);
        EXPECT_EQ(27, mixed_fe_a.get<3>()._a);
    }
//    {
//        MixedSpaceDim<Mapping> mixed_fe_a( 3, vec);
//        EXPECT_EQ(6, mixed_fe_a.get<0>()._a);
//        EXPECT_EQ(7, mixed_fe_a.get<1>()._a);
//        EXPECT_EQ(8, mixed_fe_a.get<2>()._a);
//        EXPECT_EQ(9, mixed_fe_a.get<3>()._a);
//    }
//    {
//        MixedSpaceDim<Mapping> mixed_fe_a( 3);
//        EXPECT_EQ(0, mixed_fe_a.get<0>()._a);
//        EXPECT_EQ(9, mixed_fe_a.get<1>()._a);
//        EXPECT_EQ(18, mixed_fe_a.get<2>()._a);
//        EXPECT_EQ(27, mixed_fe_a.get<3>()._a);
//    }

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
        EXPECT_EQ(13, mixed_fe.get<0>()->_a);
        EXPECT_EQ(14, mixed_fe.get<1>()->_a);
        EXPECT_EQ(15, mixed_fe.get<2>()->_a);
        EXPECT_EQ(16, mixed_fe.get<3>()->_a);
    }


    // dim and spacedim templates
    {
        MixedPtr<Mapping> mixed_fe( 3, vec);
        EXPECT_EQ(6, mixed_fe.get<0>()->_a);
        EXPECT_EQ(7, mixed_fe.get<1>()->_a);
        EXPECT_EQ(8, mixed_fe.get<2>()->_a);
        EXPECT_EQ(9, mixed_fe.get<3>()->_a);
    }
//    {
//        MixedSpaceDimPtr<Mapping> mixed_fe( 3, vec);
//        EXPECT_EQ(6, mixed_fe.get<0>()->_a);
//        EXPECT_EQ(7, mixed_fe.get<1>()->_a);
//        EXPECT_EQ(8, mixed_fe.get<2>()->_a);
//        EXPECT_EQ(9, mixed_fe.get<3>()->_a);
//    }

    // assign to base
    {
        MixedPtr<FE_XY>  fe_xy( 3, vec);
//        foo(fe_xy);
//        foo(3);
        auto mixed_fe = MixedPtr<FE>(fe_xy);
        EXPECT_EQ(13, mixed_fe.get<0>()->_a);
        EXPECT_EQ(14, mixed_fe.get<1>()->_a);
        EXPECT_EQ(15, mixed_fe.get<2>()->_a);
        EXPECT_EQ(16, mixed_fe.get<3>()->_a);
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
