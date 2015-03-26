/*
 * time_marks_test.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: jb
 */

#define DEBUG_ASSERTS_WITHOUT_MPI

#include <flow_gtest.hh>

#include "tools/time_marks.hh"


/**
 * Test for class TimeMark.
 */
TEST (TimeMark, time_mark)
{
    // Constructor.
    TimeMark tm1(-1.0, 0x01);
    TimeMark tm2(0.0, TimeMark::every_type);
    TimeMark tm3(10.0,0x05);
    TimeMark::Type my_type = 0x0a;

    //checking time values
    EXPECT_EQ(tm1.time(), -1.0);
    EXPECT_EQ(tm2.mark_type(), TimeMark::Type(~0x0));

    //checking type values - comparing masks
    EXPECT_TRUE(tm1.match_mask(0x01));
    EXPECT_TRUE(tm2.match_mask(~0x00));

    //adding type and checking mask
    tm3.add_to_type(my_type);
    EXPECT_TRUE(tm3.match_mask(0x0f));      //0x05 + 0x0a = 0x0f

    //comparing times
    EXPECT_TRUE( (tm1<tm2) && (tm2<tm3) );
}

TEST(TimeMarks, add_time_marks) {
    auto tm = TimeMarks();

    {
    auto mark_type = tm.new_mark_type();
    tm.add_time_marks(0.0, 0.1, 1.0E3, mark_type);
    auto mark_it = tm.last(mark_type);
    EXPECT_FLOAT_EQ(1.0E3, mark_it->time());
    }

    {
    auto mark_type = tm.new_mark_type();
    tm.add_time_marks(1.0, 0.1, 2.05, mark_type);

    EXPECT_FLOAT_EQ(1.0, (++tm.begin(mark_type))->time());
    EXPECT_FLOAT_EQ(2.0, tm.last(mark_type)->time());
    }
}
