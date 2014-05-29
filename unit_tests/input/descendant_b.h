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

	DescendantB();

private:
	static const int reg;
};


template <int spacedim>
const int DescendantB<spacedim>::reg =
        Input::Factory<FactoryBaseType>::template register_class< DescendantB<spacedim> >("DescendantB");


template <int spacedim>
DescendantB<spacedim>::DescendantB()
: Base<spacedim>()
{
	this->infotext_ = "Constructor of DescendantB class with spacedim = " + boost::lexical_cast<std::string>(spacedim);
};


#endif /* DESCENDANT_B_HH_ */
