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
#include <memory>
#include "type_base.hh"


namespace Input {


/**
 * @brief The Singleton class TypeRepository serves for handling the lazy-evaluated input types, derived from the base class
 * Type::TypeBase.
 *
 * When all static variables are initialized, the method TypeRepository::instance().finish() can be called
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

	/// Type stored objects of input types.
    typedef std::map< Type::TypeBase::TypeHash, std::shared_ptr<T> > TypeRepositoryMap;

    /// Public typedef of constant iterator into map of stored type.
    typedef typename TypeRepositoryMap::const_iterator TypeRepositoryMapIter;

    /// Return singleton instance of class.
    static TypeRepository & get_instance() {
    	static TypeRepository instance;
    	return instance;
    };

    /// Add @p type to TypeRepository if doesn't exist there or get existing type with same TypeHash
    std::shared_ptr<T> add_type(const T & type);

    std::shared_ptr<T> find_hash(Type::TypeBase::TypeHash hash);

    /**
     * @brief Finish all stored types.
     *
     * Iterate through all types stored in TypeRepository
     * and call finish with given status.
     *
     * Note: This method is meant to be used only with
     * FinishType::delete.
     */
    void finish(Type::FinishStatus finish_type);

    /**
     * @brief Reset and remove types marked as deleted during finish.
     *
     * Iterate through all types stored in TypeRepository and -
     *  - check count of usage of shared pointer to type (must be one)
     *  - reset this shared pointer
     *  - remove type from repository
     */
    void reset_deleted_types();

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

    /// Stores input type objects.
    TypeRepositoryMap type_repository_map_;
};


template <class T>
std::shared_ptr<T> TypeRepository<T>::add_type(const T & type) {
    Type::TypeBase::TypeHash hash = type.content_hash();

	auto search = find_hash(hash);
	if (search) {
		return search;
	} else {
		auto type_ptr = std::make_shared<T>( type );
		type_repository_map_.insert( std::pair<Type::TypeBase::TypeHash, std::shared_ptr<T>>(hash,type_ptr) );
		return type_ptr;
	}
}


template <class T>
std::shared_ptr<T> TypeRepository<T>::find_hash(Type::TypeBase::TypeHash hash) {
    auto search = type_repository_map_.find(hash);
    if (search != type_repository_map_.end()) {
        return search->second;
    } else {
        return std::shared_ptr<T>();
    }
}

template <class T>
void TypeRepository<T>::finish(Type::FinishStatus finish_type) {
	for (typename TypeRepositoryMap::iterator it = type_repository_map_.begin(); it != type_repository_map_.end(); ++it) {
		it->second->finish(finish_type);
	}
}

template <class T>
void TypeRepository<T>::reset_deleted_types() {
	std::vector< Type::TypeBase::TypeHash > deleted_hashes;
	for (typename TypeRepositoryMap::iterator it = type_repository_map_.begin(); it != type_repository_map_.end(); ++it) {
		if (it->second->finish_status() == Type::FinishStatus::delete_) {
			ASSERT(it->second.use_count() == 1)(it->second.use_count()).error();
			it->second.reset();
			deleted_hashes.push_back(it->first);
		}
	}

	for (auto deleted_hash : deleted_hashes) {
		type_repository_map_.erase(deleted_hash);
	}
}

} // namespace Input


#endif /* TYPE_REPOSITORY_HH_ */
