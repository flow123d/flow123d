/*
 * type_tuple_test.cpp
 *
 *  Created on: May 4, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>

#include <input/input_type.hh>


TEST(InputTypeTuple, declare) {
	using namespace Input::Type;

	{
		Tuple tpl = Tuple("SomeTuple", "desc.")
		   	.declare_key("start_time", Double(0.0), Default::obligatory(), "desc.")
		   	.declare_key("end_time", Double(0.0), Default::obligatory(), "desc.")
		   	.declare_key("time_step", Double(0.0), Default("0.1"), "desc.")
		   	.declare_key("time_mark", String(), Default::optional(), "desc.")
	      	.declare_key("compression", Bool(), Default("true"), "desc.")
			.declare_key("print_data", Bool(), Default("false"), "desc.")
		    .declare_key("data_description", String(), Default::optional(),"")
			.close();

		tpl.finish();

	    EXPECT_TRUE(tpl.is_finished());
	    EXPECT_EQ(tpl.auto_conversion_key_iter(), tpl.end());
	}

	{
		Tuple tpl = Tuple("OtherTuple", "desc.")
		   	.declare_key("start_time", Double(0.0), Default::obligatory(), "desc.")
		   	.declare_key("end_time", Double(0.0), Default("1.0"), "desc.")
		   	.declare_key("time_step", Double(0.0), Default("0.1"), "desc.")
		   	.declare_key("time_mark", String(), Default::optional(), "desc.")
	      	.declare_key("compression", Bool(), Default("true"), "desc.")
			.declare_key("print_data", Bool(), Default("false"), "desc.")
		    .declare_key("data_description", String(), Default::optional(),"")
			.close();

	    tpl.finish();

	    EXPECT_TRUE(tpl.is_finished());
	    EXPECT_EQ(0, (int)tpl.auto_conversion_key_iter()->key_index);
	    EXPECT_STREQ("start_time", tpl.auto_conversion_key_iter()->key_.c_str());
	}

	{
		Tuple tpl = Tuple("InvalidTuple", "desc.")
		   	.declare_key("start_time", Double(0.0), Default::obligatory(), "desc.")
		   	.declare_key("time_step", Double(0.0), Default("0.1"), "desc.")
		   	.declare_key("end_time", Double(0.0), Default::obligatory(), "desc.")
		   	.declare_key("time_mark", String(), Default::optional(), "desc.")
	      	.declare_key("compression", Bool(), Default("true"), "desc.")
		    .declare_key("data_description", String(), Default::optional(),"")
			.close();

	    EXPECT_THROW_WHAT( { tpl.finish(); }, Tuple::ExcTupleWrongKeysOrder, "in Tuple: 'InvalidTuple'" );
	}

}
