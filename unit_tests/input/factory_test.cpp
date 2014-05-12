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

#include "factory_base.h"
#include "factory_derived_a.h"
#include "factory_derived_b.h"


template class FactoryBase<3>;
template class FactoryDerivedA<3>;
template class FactoryDerivedB<3>;


TEST(FactoryTest, ClassHierarchy) {
	Input::Factory< FactoryBase<3> >::instance()->create("FactoryDerivedA", 2, 0.5);
	Input::Factory< FactoryBase<3> >::instance()->create("FactoryDerivedB");
}
