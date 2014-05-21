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

template<class Type, class... Arguments>
Factory<Type, Arguments...> * Factory<Type, Arguments...>::instance()
{
    static Factory<Type, Arguments...> factory;
    return &factory;
}


template<class Type, class... Arguments>
int Factory<Type, Arguments...>::register_function(string class_name, std::shared_ptr<Type>(* func)(Arguments...) ) {
	auto func_wrapper = std::function<std::shared_ptr<Type>(Arguments...)>(func);
	Factory<Type, Arguments...>::instance()->factory_registry_[class_name] = func_wrapper;
	return 0;
}


template<class Type, class... Arguments>
shared_ptr<Type> Factory<Type, Arguments...>::create(string name, Arguments... arguments)
{
    // find name in the registry and call factory method.
    auto it = factory_registry_.find(name);
    if(it != factory_registry_.end()) {
    	auto factory_function = boost::any_cast<std::function< shared_ptr<Type>(Arguments...)> > (it->second);
        return factory_function(arguments...);
    } else {
    	xprintf(UsrErr, "Key %s isn't registered in factory!\n", name.c_str());
    	return nullptr;
    }
}


} // namespace Input

#endif /* FACTORY_IMPL_HH_ */
