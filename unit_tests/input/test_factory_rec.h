/*
 * test_factory_rec.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef FACTORY_REC_CLASS_HH_
#define FACTORY_REC_CLASS_HH_

#include <input/type_record.hh>

#include "it_record_factory.h"


namespace it = Input::Type;

class TestFactoryRec
{
public:
	TestFactoryRec()
	{
		static it::Record i_rec = this->get_input_type().get()
				->declare_key("problem", *(ITRecordFactory::instance()->get_abstract_record("Problem").get()),
						it::Default::obligatory(), "Simulation problem to be solved.")
				.declare_key("a_tol", it::Double(0.0),
						it::Default("1.0e-9"), "Absolute residual tolerance.");
	}

	std::shared_ptr<it::Record> get_input_type()
	{
		auto type = ITRecordFactory::instance()->get_record("BaseRecord");

		return type;
	}

};

static ITRecordRegistrar rec_registrar("BaseRecord", "Base record");

#endif /* FACTORY_REC_CLASS_HH_ */
