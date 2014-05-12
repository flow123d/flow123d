/*
 * factory_derived_b.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef FACTORY_DERIVED_B_HH_
#define FACTORY_DERIVED_B_HH_

#include <boost/lexical_cast.hpp>

#include "input/factory.hh"
#include "factory_base.h"


template<int spacedim>
class FactoryDerivedB : public FactoryBase<spacedim>
{
public:
	typedef FactoryBase<spacedim> FactoryBaseType;
	static const Input::Registrar<FactoryDerivedB> reg;

	FactoryDerivedB();

};


template <int spacedim>
const Input::Registrar< FactoryDerivedB<spacedim> >
	FactoryDerivedB<spacedim>::reg("FactoryDerivedB");


template <int spacedim>
FactoryDerivedB<spacedim>::FactoryDerivedB()
: FactoryBase<spacedim>()
{
	cout << "Constructor of FactoryDerivedB class with spacedim = " << spacedim << endl;
};


#endif /* FACTORY_DERIVED_B_HH_ */
