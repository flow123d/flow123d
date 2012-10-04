/*
 * discrete_function_test.cpp
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */


/*
 * python_function_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */

#define DEBUG

#include <gtest/gtest.h>
#include <sstream>
#include <string>


#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"

#include "functions/discrete_function.hh"

using namespace std;

string input = R"CODE(
{
   discrete_1 = { constant = 1 },
   discrete_2 = 2,
   discrete_3 = { table = [1,2,3] }
}
)CODE";



TEST(DiscreteFunction, all) {
    Input::Type::Record main_rec("MainRec","");
    main_rec.declare_key("discrete_1", DiscreteFunction::get_input_type(), Input::Type::Default::obligatory(),  "");
    main_rec.declare_key("discrete_2", DiscreteFunction::get_input_type(), Input::Type::Default::obligatory(), "");
    main_rec.declare_key("discrete_3", DiscreteFunction::get_input_type(), Input::Type::Default::obligatory(), "");
    main_rec.finish();

    // read input string
    std::stringstream ss(input);
    Input::JSONToStorage reader;
    reader.read_stream( ss, main_rec);
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    DiscreteFunction fce_1( in_rec.val<Input::Record>("discrete_1") );
    EXPECT_EQ(1, fce_1.value(0));
    EXPECT_EQ(1, fce_1.value(100));

    DiscreteFunction fce_2( in_rec.val<Input::Record>("discrete_2") );
    EXPECT_EQ(2, fce_2.value(0));
    EXPECT_EQ(2, fce_2.value(100));

    DiscreteFunction fce_3( in_rec.val<Input::Record>("discrete_3") );
    EXPECT_EQ(1, fce_3.value(0));
    EXPECT_EQ(2, fce_3.value(1));
    EXPECT_EQ(3, fce_3.value(2));


    vector<unsigned int> values(3);
    values[0]=5; values[1]=6; values[2]=7;
    DiscreteFunction fce_4(values);
    EXPECT_EQ(5, fce_4.value(0));
    EXPECT_EQ(6, fce_4.value(1));
    EXPECT_EQ(7, fce_4.value(2));

}





