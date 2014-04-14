/*
 * test_factory_rec_a.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef FACTORY_REC_B_CLASS_HH_
#define FACTORY_REC_B_CLASS_HH_

#include <input/type_record.hh>

#include "it_record_factory.h"


class TestFactoryRecB
{
public:
	TestFactoryRecB()
	{
		static it::Record i_rec = this->get_input_type().get()
				->derive_from( *(ITRecordFactory::instance()->get_abstract_record("Problem").get()) )
				.declare_key("pressure_p0", it::String(),
	                    "Output stream for P0 approximation of the pressure field.");
	}

	std::shared_ptr<Input::Type::Record> get_input_type()
	{
		auto type = ITRecordFactory::instance()->get_record("UnsteadyProblem");

		return type;
	}

};

static ITRecordRegistrar unsteady_registrar("UnsteadyProblem", "Unsteady problem to solve.");

#endif /* FACTORY_REC_B_CLASS_HH_ */
