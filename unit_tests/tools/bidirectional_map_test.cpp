/*
 * bidirectional_map.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include "tools/bidirectional_map.hh"


TEST(BidirectionalMap, full_test) {
	BidirectionalMap<int> test_map;
	test_map.reserve(4);

    EXPECT_EQ(0, test_map.size());
    test_map.add_item(4);
    test_map.set_item(1, 0);
    EXPECT_EQ(1, test_map.size());
	test_map.add_item(4);
	test_map.set_item(2, 1);
	test_map.add_item(5);
	test_map.set_item(3, 2);

	EXPECT_EQ(3, test_map.size());
	for (unsigned int i=0; i<test_map.size(); ++i) {
		EXPECT_EQ(i, test_map.get_position(i+1));
		EXPECT_EQ(i+1, test_map[i]);
	}

	test_map.add_item(-6);
	EXPECT_EQ(4, test_map.size());
	EXPECT_EQ(3, test_map.get_position(-6));
	EXPECT_EQ(-6, test_map[3]);
}

