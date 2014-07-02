/*
 * factory_derived_a.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef FACTORY_DESCENDANT_A_HH_
#define FACTORY_DESCENDANT_A_HH_

#include <boost/lexical_cast.hpp>

#include "input/factory.hh"
#include "factory_base.h"

using namespace std;


template<int spacedim>
class DescendantA : public Base<spacedim>
{
public:
	typedef Base<spacedim> FactoryBaseType;

	DescendantA(int n_comp, double time);

private:
	static const int reg;
};


template <int spacedim>
const int DescendantA<spacedim>::reg =
		Input::register_class< DescendantA<spacedim>, int, double >("DescendantA");



/*
class XBase {

};

class X :public XBase {
public:
    static const int reg;
    typedef XBase FactoryBaseType;
    X(int n_comp, double time) {};
};


const int X::reg =
        Input::Factory<FactoryBaseType, int, double>::register_class< X >("DescendantA");
*/

template <int spacedim>
DescendantA<spacedim>::DescendantA(int n_comp, double time)
: Base<spacedim>()
{
	this->infotext_ = "Constructor of DescendantA class with spacedim = "
			+ boost::lexical_cast<std::string>(spacedim)
			+ ", n_comp = " + boost::lexical_cast<std::string>(n_comp)
			+ ", time = " + boost::lexical_cast<std::string>(time);
};


#endif /* FACTORY_DESCENDANT_A_HH_ */
