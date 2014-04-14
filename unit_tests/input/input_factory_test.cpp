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

#include "it_record_factory.h"
#include "test_factory_rec.h"
#include "test_factory_arec.h"
#include "test_factory_rec_a.h"
#include "test_factory_rec_b.h"



TEST(InputFactoryTest, RecordHierarchy) {
	TestFactoryRec f1;
	TestFactoryRec f2;
}
