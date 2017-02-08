/*
 * type_selection_test.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *
 *
 */


#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include <input/type_selection.hh>

using namespace Input::Type;


// Test of correct includes in type_selection.hh
TEST(InputTypeSelection, include) {
	Selection sel = Selection("SelX")
		.add_value(0, "zero")
		.close();
	EXPECT_EQ( sel.class_name(), "Selection");

}


#include <input/input_type.hh>


/**
 * Test Selection class.
 */
enum Colors {
    blue,
    white=300,
    black=45,
    red,
    green,
    yellow
};

TEST(InputTypeSelection, construction) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Selection *sel1= new Selection("Colors");
    Selection sel2=*sel1;

    sel1->add_value(blue, "blue");
    sel2.add_value(white,"white","White color");
    sel1->add_value(black,"black");
    sel2.add_value(red,"red");
    EXPECT_THROW_WHAT( {sel1->add_value(green,"red");}, feal::Exc_assert, "already exists in Selection");
    EXPECT_THROW_WHAT( {sel2.add_value(green,"red");}, feal::Exc_assert, "already exists in Selection");
    EXPECT_THROW_WHAT( {sel1->add_value(blue,"blue1");}, feal::Exc_assert, "conflicts with value");
    EXPECT_THROW_WHAT( {sel2.add_value(blue,"blue1");}, feal::Exc_assert, "conflicts with value");


    sel2.add_value(green,"green");
    sel2.close();
    EXPECT_THROW_WHAT( {sel2.add_value(yellow,"y");}, feal::Exc_assert, "in closed Selection.");

    Selection sel3;
    EXPECT_TRUE( sel3.is_finished());
    EXPECT_EQ("EmptySelection", sel3.type_name());
    EXPECT_THROW_WHAT( {sel3.add_value(1,"one");}, feal::Exc_assert, "in closed Selection.");
    // getter methods
    EXPECT_TRUE( sel2.has_name("blue") );
    EXPECT_FALSE( sel2.has_name("xblue") );
    EXPECT_TRUE( sel2.has_value(blue) );
    EXPECT_TRUE( sel2.has_value(black) );
    EXPECT_FALSE( sel2.has_value(yellow) );

    EXPECT_EQ( 45, sel2.name_to_int("black") );
    EXPECT_THROW( {sel2.name_to_int("xblack");}, Selection::ExcSelectionKeyNotFound );


    // Integer selection
    Selection int_sel = Selection("Integer selection")
    	.add_value(10, "ten")
    	.add_value(0,"zero")
    	.close();
    std::shared_ptr<Selection> sel_ptr = std::make_shared<Selection>(int_sel);


    // Selection defaults
    EXPECT_TRUE( Default("\"ten\"").check_validity( sel_ptr ) );
    EXPECT_TRUE( Default("\"zero\"").check_validity( sel_ptr ) );
    EXPECT_THROW_WHAT( { Default("\"two\"").check_validity( sel_ptr ); }, ExcWrongDefault,
            "Default value .* do not match type: " );
    EXPECT_THROW_WHAT( { Default("10").check_validity( sel_ptr ); }, ExcWrongDefault,
            "Default value .* do not match type: " );

    EXPECT_EQ(10, int_sel.name_to_int("ten"));

    delete sel1;
}
