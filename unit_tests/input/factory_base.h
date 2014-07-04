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
class Base
{
public:
	Base();

	std::string get_infotext() {
		return infotext_;
	}

protected:
	std::string infotext_;
};


template <int spacedim>
Base<spacedim>::Base()
{};


#endif /* FACTORY_BASE_HH_ */
