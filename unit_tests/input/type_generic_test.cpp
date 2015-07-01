/*
 * input_generic_test.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *
 *
 */


#include <flow_gtest.hh>

#include <input/type_generic.hh>
#include <input/type_base.hh>
#include <input/type_record.hh>


using namespace Input::Type;

TEST(GenericType, declare_key) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

	static Array parametrized_arr = Array(Parameter("parametric"));

	std::vector<TypeBase::ParameterPair> vec;
	vec.push_back( std::make_pair("parametric", boost::make_shared<Double>()) );

	boost::shared_ptr<TypeBase> type_arr = parametrized_arr.make_instance(vec);
}


static const Instance & get_generic_record(const Selection *sel, int max_limit) {
	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param1", boost::make_shared<Selection>(*sel)) );
	param_vec.push_back( std::make_pair("param2", boost::make_shared<Integer>(0, max_limit)) );

	return Instance(Record("generic_rec", "desc.")
						.declare_key("start_time", Double(), "desc.")
						.declare_key("name", String(), "desc.")
						.declare_key("param1", Parameter("param1"), "desc.")
						.declare_key("param2", Parameter("param2"), "desc.")
						.close(),
					param_vec)
			.close();
}

TEST(GenericType, generic_record) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";
	static Selection sel1 = Selection("colors")
			.add_value(0, "red")
			.add_value(1, "green")
			.add_value(2, "blue")
			.close();

	static Selection sel2 = Selection("shapes")
			.add_value(0, "rectangle")
			.add_value(1, "circle")
			.add_value(2, "triangle")
			.close();

	static Record problem = Record("problem", "desc.")
			.declare_key("primary", get_generic_record(&sel1, 10), "Primary problem.")
			.declare_key("secondary", get_generic_record(&sel2, 1000), "Secondary problem.")
			.declare_key("bool", Bool(), "Some bool key.")
			.close();

	problem.finish();
}
