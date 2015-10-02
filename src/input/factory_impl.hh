/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    factory_impl.hh
 * @brief   
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
template<class Child>
int Factory<Type, Arguments...>::register_class(string class_name) {
    auto creating_function =
        [](Arguments... args)->std::shared_ptr<Type>
        { return std::make_shared<Child>(args...); };

	auto func_wrapper = std::function<std::shared_ptr<Type>(Arguments...)>(creating_function);
	Factory<Type, Arguments...>::instance()->factory_registry_[class_name] = func_wrapper;
	//DBGMSG("Adding name: %s of class %s to factory of class %s.\n", class_name.c_str(), typeid(Child).name(), typeid(Type).name());
    return 0;
}


template<class Type, class... Arguments>
const shared_ptr<Type> Factory<Type, Arguments...>::create(string name, Arguments... arguments) const
{
    // find name in the registry and call factory method.
    auto it = factory_registry_.find(name);
    if(it != factory_registry_.end()) {
    	auto factory_function = boost::any_cast<std::function< shared_ptr<Type>(Arguments...)> > (it->second);
        return factory_function(arguments...);
    } else {
    	THROW( ExcNotRegistredClass() << EI_KeyName(name) << EI_TypeName(typeid(Type).name()) );
    	return nullptr;
    }
}


template <class ChildType, class... Arguments>
int register_class(string class_name)
{
	return Input::Factory<typename ChildType::FactoryBaseType, Arguments...>::template register_class< ChildType >(class_name);
}



} // namespace Input

#endif /* FACTORY_IMPL_HH_ */
