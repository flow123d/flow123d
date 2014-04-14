/*
 * test_factory_arec.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef FACTORY_AREC_CLASS_HH_
#define FACTORY_AREC_CLASS_HH_

#include <input/type_record.hh>

#include "it_record_factory.h"


namespace it = Input::Type;

class TestFactoryARec
{
public:
	TestFactoryARec()
	{
		cout << " call: TestFactoryARec()" << endl;
	}

	std::shared_ptr<it::AbstractRecord> get_input_type()
	{
		auto type = ITRecordFactory::instance()->get_abstract_record("Problem");

		return type;
	}

};

static ITAbstractRecordRegistrar arec_registrar("Problem", "Problem to solve.");

#endif /* FACTORY_AREC_CLASS_HH_ */
