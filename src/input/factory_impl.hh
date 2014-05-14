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
Factory<Type> * Factory<Type>::instance()
{
    static Factory<Type> factory;
    return &factory;
}


template<class Type>
template <class... Arguments>
int Factory<Type>::register_function(string class_name, std::shared_ptr<Type>(* func)(Arguments...) ) {
	auto func_wrapper = std::function<std::shared_ptr<Type>(Arguments...)>(func);
	Factory<Type>::instance()->factory_registry_[class_name] = func_wrapper;
	return 0;
}


template<class Type>
template<class... Arguments>
shared_ptr<Type> Factory<Type>::create(string name, Arguments... arguments)
{
    // find name in the registry and call factory method.
    auto it = factory_registry_.find(name);
    if(it != factory_registry_.end()) {
    	auto factory_function = boost::any_cast<std::function< shared_ptr<Type>(Arguments...)> > (it->second);
        return factory_function(arguments...);
    } else {
    	return nullptr;
    }
}


} // namespace Input

#endif /* FACTORY_IMPL_HH_ */
