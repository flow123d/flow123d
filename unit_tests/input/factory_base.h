/*
 * factory_base.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef FACTORY_BASE_HH_
#define FACTORY_BASE_HH_

#include "input/factory.hh"


template<int spacedim>
class FactoryBase
{
public:
	FactoryBase();

};

template <int spacedim>
FactoryBase<spacedim>::FactoryBase()
{};

template class Input::Factory< FactoryBase<3> >;

#endif /* FACTORY_BASE_HH_ */
