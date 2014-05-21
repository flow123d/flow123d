/*
 * input_factory_test.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */


#include <flow_gtest.hh>

#include <input/input_type.hh>
#include <input/type_record.hh>
#include <input/type_base.hh>
#include <input/type_output.hh>
#include <input/factory.hh>

#include "base.h"
#include "descendant_a.h"
#include "descendant_b.h"


template class Base<3>;
template class DescendantA<3>;
template class DescendantB<3>;


TEST(FactoryTest, ClassHierarchy) {
	EXPECT_STREQ("Constructor of DescendantA class with spacedim = 3, n_comp = 2, time = 0.5",
			( Input::Factory< Base<3>, int, double >::instance()->create("DescendantA", 2, 0.5) )->get_infotext().c_str());
	EXPECT_STREQ("Constructor of DescendantB class with spacedim = 3",
			( Input::Factory< Base<3> >::instance()->create("DescendantB") )->get_infotext().c_str());

	EXPECT_THROW_WHAT(
			{ ( Input::Factory< Base<3>, double, double >::instance() )->create("DescendantA", 1.5, 0.5); },
			ExcXprintfMsg,
			"isn't registered in factory"
			);
}
