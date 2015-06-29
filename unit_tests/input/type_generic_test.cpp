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


TEST(GenericType, declare_key) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

	static Record parametrized_rec = Record("Record","")
       	.declare_key("mesh", String(), Default("input.msh"), "Comp. mesh.")
       	.declare_key("parametric", Parameter("parametric"), Default::obligatory(), "")
		.close();

	std::vector<TypeBase::ParameterPair> vec;
	vec.push_back( std::make_pair("parametric", boost::make_shared<Double>()) );
	Instance inst = Instance(parametrized_rec, vec);

	const TypeBase &type = parametrized_rec.make_instance(vec);
}
