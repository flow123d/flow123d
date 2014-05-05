/*
 * factory_impl.hh
 *
 *  Created on: Aug 1, 2012
 *      Author: jb
 */

#ifndef FACTORY_IMPL_HH_
#define FACTORY_IMPL_HH_


namespace Input {

using namespace std;

template<class Type>
Registrar<Type>::Registrar(string name)
: factory_ref( *( Factory<BaseType>::instance() ) )
{
    // register the class factory function
	Factory<BaseType>::instance()->register_function(name, []() -> BaseType * { return new Type(); });
}


template<class Type>
Factory<Type> * Factory<Type>::instance()
{
    static Factory<Type> factory;
    return &factory;
}


template<class Type>
void Factory<Type>::register_function(string name, function<Type*(void)> class_factory_function)
{
    // register the class factory function
	field_factory_registry_[name] = class_factory_function;
}


template<class Type>
shared_ptr<Type> Factory<Type>::create(string name)
{
	Type * instance = nullptr;

    // find name in the registry and call factory method.
    auto it = field_factory_registry_.find(name);
    if(it != field_factory_registry_.end()) {
    	instance = it->second();
    }

    // wrap instance in a shared ptr and return
    if(instance != nullptr) {
    	return shared_ptr<Type>(instance);
    } else {
        return nullptr;
    }
}

} // namespace Input

#endif /* FACTORY_IMPL_HH_ */
