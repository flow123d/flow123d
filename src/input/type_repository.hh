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
 * @file    type_repository.hh
 * @brief   
 */

#ifndef TYPE_REPOSITORY_HH_
#define TYPE_REPOSITORY_HH_


#include <map>
#include "type_base.hh"


namespace Input {


/**
 * The Singleton class TypeRepository serves for handling the lazy-evaluated input types, derived from the base class
 * LazyType. When all static variables are initialized, the method TypeRepository::instance().finish() can be called
 * in order to finish initialization of lazy types such as Records, Abstracts, Arrays and Selections.
 * Selections have to be finished after all other types since they are used by Abstracts to register all
 * derived types. For this reason TypeRepository contains two arrays - one for Selections, one for the rest.
 *
 * This is list of unique instances that may contain raw pointers to possibly not yet constructed
 * (static) objects. Unique instance is the instance that creates unique instance of the data class in pimpl idiom.
 * These has to be completed/finished before use.
 *
 */
template <class T>
class TypeRepository {
public:
	/// Template parameter can be only descendant of TypeBase class.
	static_assert(std::is_base_of<Type::TypeBase, T>::value,
	        "T must be a descendant of Input::Type::TypeBase"
	    );

    typedef std::map< Type::TypeBase::TypeHash, boost::shared_ptr<T> > TypeRepositoryMap;

    /// Public typedef of constant iterator into map of stored type.
    typedef typename TypeRepositoryMap::const_iterator TypeRepositoryMapIter;

    static TypeRepository & get_instance() {
    	static TypeRepository instance;
    	return instance;
    };

    /**
     * Add @p type to TypeRepository if doesn't exist there
     * or get existing type with same TypeHash
     */
    boost::shared_ptr<T> add_type(const T & type);

    /**
     * Iterate through all types stored in TypeRepository
     * and call finish if flag @p root_of_generic_subtree_
     * has same value as @param is_root_of_generic_subtree.
     *
     * We need call finish in two steps for correct
     * functionality of generic types. In first step all
     * types marked as "root_of_generic_subtree" are
     * finished and then in second step can be finished
     * other types.
     */
    void finish(bool is_root_of_generic_subtree = false);

    /// Container-like access to the data stored in TypeRepository. Returns iterator to the first data.
    TypeRepositoryMapIter begin() const {
        return type_repository_map_.begin();
    }

    /// Container-like access to the data stored in TypeRepository. Returns iterator to the last data.
    TypeRepositoryMapIter end() const {
        return type_repository_map_.end();
    }
private:
    /// Default constructor.
    TypeRepository() {};

    TypeRepositoryMap type_repository_map_;
};


template <class T>
boost::shared_ptr<T> TypeRepository<T>::add_type(const T & type) {
    Type::TypeBase::TypeHash hash = type.content_hash();

	auto search = type_repository_map_.find(hash);
	if (search != type_repository_map_.end()) {
		return search->second;
	} else {
		auto type_ptr = boost::make_shared<T>( type );
		type_repository_map_.insert( std::pair<Type::TypeBase::TypeHash, boost::shared_ptr<T>>(hash,type_ptr) );
		return type_ptr;
	}
}

template <class T>
void TypeRepository<T>::finish(bool is_root_of_generic_subtree) {
	// We need reverse iterating for correct finish of generic types.
	for (typename TypeRepositoryMap::reverse_iterator it = type_repository_map_.rbegin(); it != type_repository_map_.rend(); ++it) {
		if (is_root_of_generic_subtree == it->second->is_root_of_generic_subtree()) {
			if (is_root_of_generic_subtree) {
#ifdef FLOW123D_DEBUG
           		if ( !it->second->is_finished() ) xprintf(Warn, "Unused root of generic subtree: '%s'.\n", it->second->type_name().c_str());
#endif
				it->second->finish(true);
			} else {
				it->second->finish();
			}
		}
	}
}

} // namespace Input


#endif /* TYPE_REPOSITORY_HH_ */
