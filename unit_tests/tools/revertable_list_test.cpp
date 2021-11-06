/*
 * tmp_size_list_test.cpp
 *
 *  Created on: Jul 27, 2020
 *      Author: df
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include "tools/revertable_list.hh"


TEST(RevertableList, all) {
    std::array<int, 10> data = {0, 2, 4, 6, 8, 1, 3, 5, 7, 9};

    RevertableList<int> list(10);
    EXPECT_EQ(list.permanent_size(), 0);
    EXPECT_EQ(list.temporary_size(), 0);
    EXPECT_EQ(list.reserved_size(), 10);

    for (uint i=0; i<5; ++i)
        list.push_back( data[i] );
    EXPECT_EQ(list.permanent_size(), 0);
    EXPECT_EQ(list.temporary_size(), 5);

    list.make_permanent();
    EXPECT_EQ(list.permanent_size(), 5);
    EXPECT_EQ(list.temporary_size(), 5);

    for (uint i=5; i<10; ++i)
        list.push_back( data[i] );
    EXPECT_EQ(list.temporary_size(), 10);

    list.revert_temporary();
    EXPECT_EQ(list.temporary_size(), 5);

    for (uint i=0; i<5; ++i)
        EXPECT_EQ(list[i], data[i]);
}
