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

using namespace std;


template<int spacedim>
class FactoryDerivedA : public FactoryBase<spacedim>
{
public:
	typedef FactoryBase<spacedim> FactoryBaseType;
	static const int reg;

	static shared_ptr< FactoryBase<spacedim> > create_instance(int n_comp, double time) {
		return make_shared< FactoryDerivedA<spacedim> >(n_comp, time);
	}

	FactoryDerivedA(int n_comp, double time);

};


template <int spacedim>
const int FactoryDerivedA<spacedim>::reg =
		Input::Factory<FactoryBaseType>::register_function("FactoryDerivedA", FactoryDerivedA<spacedim>::create_instance );


template <int spacedim>
FactoryDerivedA<spacedim>::FactoryDerivedA(int n_comp, double time)
: FactoryBase<spacedim>()
{
	cout << "Constructor of FactoryDerivedA class with spacedim = " << spacedim;
	cout << ", n_comp = " << n_comp << ", time = " << time << endl;
};


#endif /* FACTORY_DERIVED_A_HH_ */
