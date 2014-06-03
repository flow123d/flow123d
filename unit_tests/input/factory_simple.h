/*
 * factory_simple.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef FACTORY_SIMPLE_HH_
#define FACTORY_SIMPLE_HH_

#include <boost/lexical_cast.hpp>

#include "input/factory.hh"


using namespace std;


class SimpleBase
{
public:
	SimpleBase() {}

	std::string get_infotext() {
		return infotext_;
	}

protected:
	std::string infotext_;
};


class SimpleDescendant : public SimpleBase {
public:
    static const int reg;
    typedef SimpleBase FactoryBaseType;
    SimpleDescendant(int n_comp, double time);
};


const int SimpleDescendant::reg =
        Input::Factory<FactoryBaseType, int, double>::register_class< SimpleDescendant >("SimpleDescendant");

SimpleDescendant::SimpleDescendant(int n_comp, double time) {
	this->infotext_ = "Constructor of SimpleDescendant class with n_comp = "
			+ boost::lexical_cast<std::string>(n_comp) + ", time = "
			+ boost::lexical_cast<std::string>(time);
};


#endif /* FACTORY_SIMPLE_HH_ */
