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
 * @file    type_generic.cc
 * @brief   
 */

#include <input/type_generic.hh>
#include <input/type_repository.hh>
#include <input/attribute_lib.hh>
#include "system/asserts.hh"                           // for Assert, ASSERT

#include <boost/functional/hash.hpp>



namespace Input {

namespace Type {


/*******************************************************************
 * implementation of Parameter
 */

Parameter::Parameter(const string & parameter_name)
: name_(parameter_name) {}


Parameter::Parameter(const Parameter & other)
: TypeBase(other), name_(other.name_) {}


string Parameter::type_name() const {
    return name_;
}


string Parameter::class_name() const {
	return "Parameter";
}


TypeBase::TypeHash Parameter::content_hash() const {
	TypeHash seed=0;
	boost::hash_combine(seed, "Parameter");
	boost::hash_combine(seed, type_name());

	return seed;
}


TypeBase::MakeInstanceReturnType Parameter::make_instance(std::vector<ParameterPair> vec) {

    // Find the parameter value in the incoming vector.
	auto parameter_iter = std::find_if(vec.begin(), vec.end(),
	                                   [this](const ParameterPair & item) -> bool { return item.first == this->name_; });
	if (parameter_iter != vec.end()) {
	    ParameterMap parameter_map;
		parameter_map[parameter_iter->first] = parameter_iter->second->content_hash();
		return std::make_pair( parameter_iter->second, parameter_map );
	} else {
	    // throw if the parameter value is missing
	    THROW( ExcParamaterNotSubsituted() << EI_Object(this->name_));
	}
}


FinishStatus Parameter::finish(FinishStatus finish_type) {
	ASSERT(finish_type != FinishStatus::none_).error();

	if (finish_type == FinishStatus::regular_) THROW( ExcParamaterInIst() << EI_Object(this->name_));
	return finish_type;
}


/*******************************************************************
 * implementation of Instance
 */

Instance::Instance(TypeBase &generic_type, std::vector<TypeBase::ParameterPair> parameters)
: generic_type_(generic_type), parameters_(parameters) {}


TypeBase::TypeHash Instance::content_hash() const {
	TypeHash seed=0;
	boost::hash_combine(seed, "Instance");
	boost::hash_combine(seed, generic_type_.content_hash() );
	for (std::vector<TypeBase::ParameterPair>::const_iterator it = parameters_.begin(); it!=parameters_.end(); it++) {
		boost::hash_combine(seed, (*it).first );
		boost::hash_combine(seed, (*it).second->content_hash() );
	}

	return seed;
}


const Instance &Instance::close() const {
	return *( Input::TypeRepository<Instance>::get_instance().add_type( *this ) );
}


FinishStatus Instance::finish(FinishStatus finish_type) {
	return generic_type_.finish(finish_type);
}


/// Print parameter vector to formatted string.


// Implements @p TypeBase::make_instance.
TypeBase::MakeInstanceReturnType Instance::make_instance(std::vector<ParameterPair> vec) {
	// check if instance is created
	if (created_instance_.first) {
		return created_instance_;
	}

	try {
		created_instance_ = generic_type_.make_instance(parameters_);
	} catch (ExcParamaterNotSubsituted &e) {

	    ParameterMap aux_map;
	    for(auto &item : vec) aux_map[item.first]=0;
        e << EI_ParameterList( TypeBase::print_parameter_map_keys_to_json(aux_map) );
        throw;
	}




#ifdef FLOW123D_DEBUG_ASSERTS
	for (std::vector<TypeBase::ParameterPair>::const_iterator vec_it = parameters_.begin(); vec_it!=parameters_.end(); vec_it++) {
		ParameterMap::iterator map_it = created_instance_.second.find( vec_it->first );

        ParameterMap aux_map;
        for(auto &item : vec) aux_map[item.first]=0;

		ASSERT_DBG(map_it != created_instance_.second.end())(vec_it->first)(generic_type_.type_name())
				.error("Unused parameter in input type instance");
	}
#endif
	return created_instance_;
}

} // closing namespace Type
} // closing namespace Input
