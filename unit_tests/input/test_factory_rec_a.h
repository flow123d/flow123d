/*
 * test_factory_rec_a.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef FACTORY_REC_A_CLASS_HH_
#define FACTORY_REC_A_CLASS_HH_

#include <input/type_record.hh>

#include "it_record_factory.h"


class TestFactoryRecA
{
public:
	TestFactoryRecA()
	{
		static it::Record i_rec = this->get_input_type().get()
				->derive_from( *(ITRecordFactory::instance()->get_abstract_record("Problem").get()) )
				.declare_key("time", it::Double(0.0), it::Default("1.0e-9"),
						"Apply field setting in this record after this time.\n"
						"These times have to form an increasing sequence.");
	}

	std::shared_ptr<Input::Type::Record> get_input_type()
	{
		auto type = ITRecordFactory::instance()->get_record("SteadyProblem");

		return type;
	}

};

static ITRecordRegistrar steady_registrar("SteadyProblem", "Steady problem to solve.");

#endif /* FACTORY_REC_A_CLASS_HH_ */
