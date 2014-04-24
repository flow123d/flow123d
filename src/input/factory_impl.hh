/*
 * factory_impl.hh
 *
 *  Created on: Aug 1, 2012
 *      Author: jb
 */

#ifndef FACTORY_IMPL_HH_
#define FACTORY_IMPL_HH_

#include "input/factory.hh"
#include "fields/field_base.hh"


namespace Input {

using namespace std;

template<class Type>
Registrar<Type>::Registrar(string name)
{
    // register the class factory function
	Factory::instance()->register_function<Type>(name, [](void) -> Type * { return new Type();});
}


template<class Type>
void Factory::register_function(string name, function<Type*(void)> class_factory_function)
{
    // register the class factory function
	field_factory_registry_[name] = class_factory_function;
}


template<class Type>
shared_ptr<Type> Factory::create(AbstractRecord &rec)
{
	Type * instance = nullptr;
	string name = rec.type().type_name();

    // find name in the registry and call factory method.
    auto it = field_factory_registry_.find(name);
    if(it != field_factory_registry_.end())
        instance = it->second();

    // wrap instance in a shared ptr and return
    if(instance != nullptr)
        return shared_ptr<Type>(instance);
    else
        return nullptr;
}

} // namespace Input

#endif /* FACTORY_IMPL_HH_ */
