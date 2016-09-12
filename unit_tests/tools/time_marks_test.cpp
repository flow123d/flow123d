/*
 * time_marks_test.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: jb
 */

#define DEBUG_ASSERTS_WITHOUT_MPI
#define FEAL_OVERRIDE_ASSERTS

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
     TimeMarks tm;

    {
    auto mark_type = tm.new_mark_type();
    tm.add_time_marks(0.0, 0.1, 1.0E3, mark_type);
    auto mark_it = tm.last(mark_type);
    EXPECT_FLOAT_EQ(1.0E3, mark_it->time());
    EXPECT_EQ(mark_type, mark_it->mark_type());
    TimeMark changed_tm = tm.add(TimeMark(0.0, tm.type_output()));
    EXPECT_EQ(mark_type | tm.type_output(), changed_tm.mark_type());
    mark_it = tm.begin(mark_type); ++mark_it;

    EXPECT_EQ(mark_type | tm.type_output(), mark_it->mark_type());
    }

    {
    auto mark_type = tm.new_mark_type();
    tm.add_time_marks(1.0, 0.1, 2.05, mark_type);

    EXPECT_FLOAT_EQ(1.0, (++tm.begin(mark_type))->time());
    EXPECT_FLOAT_EQ(2.0, tm.last(mark_type)->time());

    auto mark_type_2 = tm.new_mark_type();
    tm.add_to_type_all(mark_type, mark_type_2);

    TimeMark time_mark = *(tm.begin(mark_type));
    EXPECT_TRUE( time_mark.match_mask(mark_type_2) );
    EXPECT_FLOAT_EQ(time_mark.time(), tm.begin(mark_type_2)->time());
    }
}
