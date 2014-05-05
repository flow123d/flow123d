/*
 * factory_derived_a.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef FACTORY_DERIVED_A_HH_
#define FACTORY_DERIVED_A_HH_

#include <boost/lexical_cast.hpp>

#include "input/factory.hh"
#include "factory_base.h"


template<int spacedim>
class FactoryDerivedA : public FactoryBase<spacedim>
{
public:
	typedef FactoryBase<spacedim> FactoryBaseType;
	static const Input::Registrar<FactoryDerivedA> reg;

	FactoryDerivedA();

};

template <int spacedim>
const Input::Registrar< FactoryDerivedA<spacedim> > FactoryDerivedA<spacedim>::reg =
		Input::Registrar< FactoryDerivedA<spacedim> >("FactoryDerivedA" + boost::lexical_cast<std::string>(spacedim));


template <int spacedim>
FactoryDerivedA<spacedim>::FactoryDerivedA()
: FactoryBase<spacedim>()
{
	cout << "Constructor of FactoryDerivedA class with spacedim = " << spacedim << endl;
};


#endif /* FACTORY_DERIVED_A_HH_ */
