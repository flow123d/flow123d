/*
 * mixed_test.cpp
 *
 *  Created on: Apr 8, 2019
 *      Author: jb
 */

#include <flow_gtest.hh>
#include <vector>
#include "tools/mixed.hh"

template <int dim, int spacedim>
class Mapping {
public:
    Mapping(int a, std::vector<int> b) : _a(a + dim + spacedim) {
    }

    Mapping(int a) : _a(a * dim * spacedim) {
    }

    int _a;
};


template <int dim>
class FE {
public:
    FE(int a, std::vector<int> b) : _a(a + dim + 10) {
    }

    FE(int a) : _a(a * dim * 10) {
    }

    int _a;
};


TEST(Mixed, mixed) {
    // only dim parameter templates
    std::vector<int> vec = {1, 2};
    {
        Mixed<FE> mixed_fe_a = Mixed<FE>( 3, vec);
        EXPECT_EQ(13, mixed_fe_a.get<0>()._a);
        EXPECT_EQ(14, mixed_fe_a.get<1>()._a);
        EXPECT_EQ(15, mixed_fe_a.get<2>()._a);
        EXPECT_EQ(16, mixed_fe_a.get<3>()._a);
    }
    {
        Mixed<FE> mixed_fe_a = Mixed<FE>( 3);
        EXPECT_EQ(0, mixed_fe_a.get<0>()._a);
        EXPECT_EQ(30, mixed_fe_a.get<1>()._a);
        EXPECT_EQ(60, mixed_fe_a.get<2>()._a);
        EXPECT_EQ(90, mixed_fe_a.get<3>()._a);
    }


    // dim and spacedim templates
    {
        Mixed<FixSpaceDim<Mapping>::template type> mixed_fe_a( 3, vec);
        EXPECT_EQ(6, mixed_fe_a.get<0>()._a);
        EXPECT_EQ(7, mixed_fe_a.get<1>()._a);
        EXPECT_EQ(8, mixed_fe_a.get<2>()._a);
        EXPECT_EQ(9, mixed_fe_a.get<3>()._a);
    }
    {
        Mixed<FixSpaceDim<Mapping>::template type> mixed_fe_a( 3);
        EXPECT_EQ(0, mixed_fe_a.get<0>()._a);
        EXPECT_EQ(9, mixed_fe_a.get<1>()._a);
        EXPECT_EQ(18, mixed_fe_a.get<2>()._a);
        EXPECT_EQ(27, mixed_fe_a.get<3>()._a);
    }
    {
        MixedSpaceDim<Mapping> mixed_fe_a( 3, vec);
        EXPECT_EQ(6, mixed_fe_a.get<0>()._a);
        EXPECT_EQ(7, mixed_fe_a.get<1>()._a);
        EXPECT_EQ(8, mixed_fe_a.get<2>()._a);
        EXPECT_EQ(9, mixed_fe_a.get<3>()._a);
    }
    {
        MixedSpaceDim<Mapping> mixed_fe_a( 3);
        EXPECT_EQ(0, mixed_fe_a.get<0>()._a);
        EXPECT_EQ(9, mixed_fe_a.get<1>()._a);
        EXPECT_EQ(18, mixed_fe_a.get<2>()._a);
        EXPECT_EQ(27, mixed_fe_a.get<3>()._a);
    }

}

