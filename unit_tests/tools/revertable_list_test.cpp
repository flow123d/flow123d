/*
 * tmp_size_list_test.cpp
 *
 *  Created on: Jul 27, 2020
 *      Author: df
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include "tools/revertable_list.hh"


TEST(RevertableList, fixed_reserved_size) {
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

TEST(RevertableList, growing_reserved_size) {
    std::array<int, 10> data = {0, 2, 4, 6, 8, 1, 3, 5, 7, 9};

    RevertableList<int> list(2, 2);
    EXPECT_EQ(list.permanent_size(), 0);
    EXPECT_EQ(list.temporary_size(), 0);
    EXPECT_EQ(list.reserved_size(), 2);

    for (uint i=0; i<5; ++i)
        list.push_back( data[i] );
    EXPECT_EQ(list.permanent_size(), 0);
    EXPECT_EQ(list.temporary_size(), 5);
    EXPECT_EQ(list.reserved_size(), 6);

    list.make_permanent();
    EXPECT_EQ(list.permanent_size(), 5);
    EXPECT_EQ(list.temporary_size(), 5);

    for (uint i=5; i<10; ++i)
        list.push_back( data[i] );
    EXPECT_EQ(list.temporary_size(), 10);
    EXPECT_EQ(list.reserved_size(), 10);

    list.revert_temporary();
    EXPECT_EQ(list.temporary_size(), 5);

    for (uint i=0; i<5; ++i)
        EXPECT_EQ(list[i], data[i]);
}

TEST(RevertableList, default_constructor) {
    std::array<int, 10> data = {0, 2, 4, 6, 8, 1, 3, 5, 7, 9};

    // base usage of default constructor
    RevertableList<int> list;
    EXPECT_EQ(list.permanent_size(), 0);
    EXPECT_EQ(list.temporary_size(), 0);
    EXPECT_EQ(list.reserved_size(), 0);

    for (uint i=0; i<5; ++i) {
        list.push_back( data[i] );
        EXPECT_EQ(list.permanent_size(), 0);
        EXPECT_EQ(list.temporary_size(), i+1);
        // capacity of list is resized automatically, test of reserved size is irrelevant
    }

    list.make_permanent();
    EXPECT_EQ(list.permanent_size(), 5);
    EXPECT_EQ(list.temporary_size(), 5);

    for (uint i=5; i<10; ++i) {
        list.push_back( data[i] );
        EXPECT_EQ(list.temporary_size(), i+1);
    }

    list.revert_temporary();
    EXPECT_EQ(list.temporary_size(), 5);

    for (uint i=0; i<5; ++i)
        EXPECT_EQ(list[i], data[i]);

    // modification by reinit_default_list method
    RevertableList<int> list2;
    EXPECT_EQ(list2.permanent_size(), 0);
    EXPECT_EQ(list2.temporary_size(), 0);
    EXPECT_EQ(list2.reserved_size(), 0);

    list2.reinit_default_list(5, 5);
    for (uint i=0; i<6; ++i)
        list2.push_back( data[i] );
    EXPECT_EQ(list2.permanent_size(), 0);
    EXPECT_EQ(list2.temporary_size(), 6);
    EXPECT_EQ(list2.reserved_size(), 10);
    list2.make_permanent();
    EXPECT_EQ(list2.permanent_size(), 6);
}
