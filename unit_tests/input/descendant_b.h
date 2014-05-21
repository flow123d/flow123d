/*
 * factory_derived_b.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef DESCENDANT_B_HH_
#define DESCENDANTD_B_HH_

#include <boost/lexical_cast.hpp>

#include "input/factory.hh"
#include "base.h"


template<int spacedim>
class DescendantB : public Base<spacedim>
{
public:
	typedef Base<spacedim> FactoryBaseType;

	static shared_ptr< Base<spacedim> > create_instance() {
		return make_shared< DescendantB<spacedim> >();
	}

	DescendantB();

private:
	static const int reg;
};


template <int spacedim>
const int DescendantB<spacedim>::reg =
		Input::Factory<FactoryBaseType>::register_function("DescendantB", DescendantB<spacedim>::create_instance );


template <int spacedim>
DescendantB<spacedim>::DescendantB()
: Base<spacedim>()
{
	this->infotext_ = "Constructor of DescendantB class with spacedim = " + boost::lexical_cast<std::string>(spacedim);
};


#endif /* DESCENDANT_B_HH_ */
