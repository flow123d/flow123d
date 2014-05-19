/*
 * base.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef BASE_HH_
#define BASE_HH_

#include "input/factory.hh"


template<int spacedim>
class Base
{
public:
	Base();

	static Input::Type::AbstractRecord input_type;

	std::string get_infotext() {
		return infotext_;
	}

protected:
	std::string infotext_;
};


template <int spacedim>
Base<spacedim>::Base()
{};


template <int spacedim>
Input::Type::AbstractRecord Base<spacedim>::input_type
    = Input::Type::AbstractRecord("Base:"+spacedim, "Abstract record desc.");


//template class Input::Factory< Base<3> >;

#endif /* BASE_HH_ */
