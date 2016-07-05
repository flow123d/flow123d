/*
 * ref_element_test.cpp
 *
 *  Created on: Jul 4, 2016
 *      Author: jb
 */



#include <flow_gtest.hh>
#include "mesh/ref_element.hh"

#include "system/armadillo_setup.hh"
#include "arma_expect.hh"
#include "system/sys_profiler.hh"
#include "armadillo"

using namespace std;

TEST(RefElement, barycentric_on_face) {
    armadillo_setup();

    // dim 1
    EXPECT_ARMA_EQ( arma::vec("1"),
            RefElement<1>::barycentric_on_face( arma::vec("0 1"), 0));
    EXPECT_ARMA_EQ( arma::vec("1"),
            RefElement<1>::barycentric_on_face( arma::vec("1 0"), 1));

    // dim 2
    EXPECT_ARMA_EQ( arma::vec("0 1"),             RefElement<2>::barycentric_on_face( arma::vec("0 0 1"), 0));
    EXPECT_ARMA_EQ( arma::vec("1 0"),             RefElement<2>::barycentric_on_face( arma::vec("1 0 0"), 0));

    EXPECT_ARMA_EQ( arma::vec("0 1"),             RefElement<2>::barycentric_on_face( arma::vec("0 0 1"), 1));
    EXPECT_ARMA_EQ( arma::vec("1 0"),             RefElement<2>::barycentric_on_face( arma::vec("0 1 0"), 1));

    EXPECT_ARMA_EQ( arma::vec("0 1"),             RefElement<2>::barycentric_on_face( arma::vec("1 0 0"), 2));
    EXPECT_ARMA_EQ( arma::vec("1 0"),             RefElement<2>::barycentric_on_face( arma::vec("0 1 0"), 2));


    // dim 3
    EXPECT_ARMA_EQ( arma::vec("0 0 1"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 0 1"), 0));
    EXPECT_ARMA_EQ( arma::vec("1 0 0"),             RefElement<3>::barycentric_on_face( arma::vec("1 0 0 0"), 0));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),             RefElement<3>::barycentric_on_face( arma::vec("0 1 0 0"), 0));

    EXPECT_ARMA_EQ( arma::vec("0 0 1"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 0 1"), 1));
    EXPECT_ARMA_EQ( arma::vec("1 0 0"),             RefElement<3>::barycentric_on_face( arma::vec("1 0 0 0"), 1));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 1 0"), 1));

    EXPECT_ARMA_EQ( arma::vec("0 0 1"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 0 1"), 2));
    EXPECT_ARMA_EQ( arma::vec("1 0 0"),             RefElement<3>::barycentric_on_face( arma::vec("0 1 0 0"), 2));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 1 0"), 2));

    EXPECT_ARMA_EQ( arma::vec("0 0 1"),             RefElement<3>::barycentric_on_face( arma::vec("1 0 0 0"), 3));
    EXPECT_ARMA_EQ( arma::vec("1 0 0"),             RefElement<3>::barycentric_on_face( arma::vec("0 1 0 0"), 3));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 1 0"), 3));
}

TEST(RefElement, barycentric_from_face) {
    armadillo_setup();

    // dim 1
    EXPECT_ARMA_EQ( arma::vec("0 1"),
            RefElement<1>::barycentric_from_face( arma::vec("1"), 0));
    EXPECT_ARMA_EQ( arma::vec("1 0"),
            RefElement<1>::barycentric_from_face( arma::vec("1"), 1));

    // dim 2
    EXPECT_ARMA_EQ( arma::vec("0 0 1"),
            RefElement<2>::barycentric_from_face( arma::vec("0 1"), 0));
    EXPECT_ARMA_EQ( arma::vec("1 0 0"),
            RefElement<2>::barycentric_from_face( arma::vec("1 0"), 0));

    EXPECT_ARMA_EQ( arma::vec("0 0 1"),
            RefElement<2>::barycentric_from_face( arma::vec("0 1"), 1));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),
            RefElement<2>::barycentric_from_face( arma::vec("1 0"), 1));

    EXPECT_ARMA_EQ( arma::vec("1 0 0"),
            RefElement<2>::barycentric_from_face( arma::vec("0 1"), 2));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),
            RefElement<2>::barycentric_from_face( arma::vec("1 0"), 2));


    // dim 3
    EXPECT_ARMA_EQ( arma::vec("0 0 0 1"),
            RefElement<3>::barycentric_from_face( arma::vec("0 0 1"), 0));
    EXPECT_ARMA_EQ( arma::vec("1 0 0 0"),
            RefElement<3>::barycentric_from_face( arma::vec("1 0 0"), 0));
    EXPECT_ARMA_EQ( arma::vec("0 1 0 0"),
            RefElement<3>::barycentric_from_face( arma::vec("0 1 0"), 0));

    EXPECT_ARMA_EQ( arma::vec("0 0 0 1"),
            RefElement<3>::barycentric_from_face( arma::vec("0 0 1"), 1));
    EXPECT_ARMA_EQ( arma::vec("1 0 0 0"),
            RefElement<3>::barycentric_from_face( arma::vec("1 0 0"), 1));
    EXPECT_ARMA_EQ( arma::vec("0 0 1 0"),
            RefElement<3>::barycentric_from_face( arma::vec("0 1 0"), 1));

    EXPECT_ARMA_EQ( arma::vec("0 0 0 1"),
            RefElement<3>::barycentric_from_face( arma::vec("0 0 1"), 2));
    EXPECT_ARMA_EQ( arma::vec("0 1 0 0"),
            RefElement<3>::barycentric_from_face( arma::vec("1 0 0"), 2));
    EXPECT_ARMA_EQ( arma::vec("0 0 1 0"),
            RefElement<3>::barycentric_from_face( arma::vec("0 1 0"), 2));

    EXPECT_ARMA_EQ( arma::vec("1 0 0 0"),
            RefElement<3>::barycentric_from_face( arma::vec("0 0 1"), 3));
    EXPECT_ARMA_EQ( arma::vec("0 1 0 0"),
            RefElement<3>::barycentric_from_face( arma::vec("1 0 0"), 3));
    EXPECT_ARMA_EQ( arma::vec("0 0 1 0"),
            RefElement<3>::barycentric_from_face( arma::vec("0 1 0"), 3));
}

TEST(RefElement, centers_of_subelements) {
    armadillo_setup();

    // dim = 1
    {
    auto list = RefElement<1>::centers_of_subelements(0);
    EXPECT_EQ( 2, list.size());
    EXPECT_ARMA_EQ( arma::vec("0"), list[0]);
    EXPECT_ARMA_EQ( arma::vec("1"), list[1]);
    }
    {
    auto list = RefElement<1>::centers_of_subelements(1);
    EXPECT_EQ( 1, list.size());
    EXPECT_ARMA_EQ( arma::vec("0.5"), list[0]);
    }

    // dim = 2
    {
    auto list = RefElement<2>::centers_of_subelements(0);
    EXPECT_EQ( 3, list.size());
    EXPECT_ARMA_EQ( arma::vec("0 0"), list[0]);
    EXPECT_ARMA_EQ( arma::vec("1 0"), list[1]);
    EXPECT_ARMA_EQ( arma::vec("0 1"), list[2]);
    }
    {
    auto list = RefElement<2>::centers_of_subelements(1);
    EXPECT_EQ( 3, list.size());
    EXPECT_ARMA_EQ( arma::vec("0.5 0"), list[0]);
    EXPECT_ARMA_EQ( arma::vec("0 0.5"), list[1]);
    EXPECT_ARMA_EQ( arma::vec("0.5 0.5"), list[2]);
    }
    {
    auto list = RefElement<2>::centers_of_subelements(2);
    EXPECT_EQ( 1, list.size());
    EXPECT_ARMA_EQ( arma::vec({ 1/3.0, 1/3.0 }), list[0]);
    }

}


TEST(RefElement, clip_1d) {
    armadillo_setup();

    // in element
    EXPECT_ARMA_EQ( arma::vec("0 1"), RefElement<1>::clip( arma::vec("0 1")));
    EXPECT_ARMA_EQ( arma::vec("0.5 0.5"), RefElement<1>::clip( arma::vec("0.5 0.5")));
    EXPECT_ARMA_EQ( arma::vec("1 0"), RefElement<1>::clip( arma::vec("1 0")));
    // out of element
    EXPECT_ARMA_EQ( arma::vec("0 1"), RefElement<1>::clip( arma::vec("-0.5 1.5")));
    EXPECT_ARMA_EQ( arma::vec("1 0"), RefElement<1>::clip( arma::vec("1.5 -0.5")));
}


TEST(RefElement, clip_2d) {
    armadillo_setup();

    //in element
    EXPECT_ARMA_EQ( arma::vec("1 0 0"), RefElement<2>::clip( arma::vec("1 0 0")));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"), RefElement<2>::clip( arma::vec("0 1 0")));
    EXPECT_ARMA_EQ( arma::vec("0 0 1"), RefElement<2>::clip( arma::vec("0 0 1")));
    EXPECT_ARMA_EQ( arma::vec("0.3 0.3 0.4"), RefElement<2>::clip( arma::vec("0.3 0.3 0.4")));
    // out of element
    EXPECT_ARMA_EQ( arma::vec("0 0 1"), RefElement<2>::clip( arma::vec("-0.5 0 1.5")));
    EXPECT_ARMA_EQ( arma::vec("0 0 1"), RefElement<2>::clip( arma::vec("0 -0.5 1.5")));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"), RefElement<2>::clip( arma::vec("-0.5 1 0.5")));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"), RefElement<2>::clip( arma::vec("0.3 1.3 -0.6")));
    EXPECT_ARMA_EQ( arma::vec("1 0 0"), RefElement<2>::clip( arma::vec("1 -0.5 0.5")));
    EXPECT_ARMA_EQ( arma::vec("1 0 0"), RefElement<2>::clip( arma::vec("1.3 0.3 -0.6")));

    EXPECT_ARMA_EQ( arma::vec("0.5 0 0.5"), RefElement<2>::clip( arma::vec("0.5 -0.3 0.8")));
    EXPECT_ARMA_EQ( arma::vec("0 0.5 0.5"), RefElement<2>::clip( arma::vec("-0.3 0.5 0.8")));
    EXPECT_ARMA_EQ( arma::vec("0.5 0.5 0"), RefElement<2>::clip( arma::vec("1 1 -1")));
}

