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
 * @file    type_base.cc
 * @brief   
 */

#include <limits>
#include <ios>
#include <map>
#include <vector>
#include <string>
#include <iomanip>

#include "system/system.hh"

#include <boost/type_traits.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include <boost/pointer_cast.hpp>


#include "input_type.hh"
#include "type_output.hh"
#include "type_repository.hh"
#include "attribute_lib.hh"
#include "json_spirit/json_spirit.h"


namespace Input {
namespace Type {

using namespace std;



/*******************************************************************
 * implementation of TypeBase
 */



TypeBase::TypeBase()
: attributes_( boost::make_shared<attribute_map>() ), root_of_generic_subtree_(false),
  generic_type_hash_(0) {}



TypeBase::TypeBase(const TypeBase& other)
: attributes_(other.attributes_), root_of_generic_subtree_(other.root_of_generic_subtree_),
  generic_type_hash_(other.generic_type_hash_), parameter_map_(other.parameter_map_) {}



TypeBase::~TypeBase() {}


bool TypeBase::is_valid_identifier(const string& key) {
  namespace ba = boost::algorithm;
  return ba::all( key, ba::is_lower() || ba::is_digit() || ba::is_any_of("_") );
}


string TypeBase::desc() const {
    stringstream ss;
    ss << OutputText(this);
    return ss.str();
}



void TypeBase::lazy_finish() {
	Input::TypeRepository<Instance>::get_instance().finish();
	Input::TypeRepository<Abstract>::get_instance().finish(true);
	Input::TypeRepository<Record>::get_instance().finish(true);
	Input::TypeRepository<Abstract>::get_instance().finish();
	Input::TypeRepository<Record>::get_instance().finish();
	Input::TypeRepository<Selection>::get_instance().finish();
}


 std::string TypeBase::hash_str(TypeHash hash) {
    stringstream ss;
    ss << "\"" << std::hex << hash << "\"";
    return ss.str();
}




void TypeBase::add_attribute(std::string name, json_string val) {
	if (validate_json(val)) {
		(*attributes_)[name] = val;
	} else {
		xprintf(PrgErr, "Invalid JSON format of attribute '%s'.\n", name.c_str());
	}
}


bool TypeBase::validate_json(json_string str) const {
    try {
    	json_spirit::mValue node;
    	json_spirit::read_or_throw( str, node);
    	return true;
    } catch (json_spirit::Error_position &e ) {
        return false;
    }
}


TypeBase::json_string TypeBase::print_parameter_map_to_json(ParameterMap parameter_map) const {
	std::stringstream ss;
	ss << "[";
	for (ParameterMap::iterator it=parameter_map.begin(); it!=parameter_map.end(); it++) {
		if (it != parameter_map.begin()) ss << "," << endl;
		ss << "{ \"" << (*it).first << "\" : " << TypeBase::hash_str( (*it).second ) << " }";
	}
	ss << "]";
	return ss.str();
}


void TypeBase::set_parameters_attribute(ParameterMap parameter_map) {
	this->add_attribute(FlowAttributes::parameters(), this->print_parameter_map_to_json(parameter_map));
}










std::ostream& operator<<(std::ostream& stream, const TypeBase& type) {
    return ( stream << OutputText(&type) );
}



/**********************************************************************************
 * implementation of Type::Array
 */

TypeBase::TypeHash Array::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, data_->lower_bound_);
    boost::hash_combine(seed, data_->upper_bound_);
    boost::hash_combine(seed, data_->type_of_values_->content_hash() );
    return seed;
}


bool Array::finish(bool is_generic) {
	return data_->finish(is_generic);
}



bool Array::ArrayData::finish(bool is_generic)
{
	if (finished) return true;

	if (typeid( *(type_of_values_.get()) ) == typeid(Instance)) type_of_values_ = type_of_values_->make_instance().first;
	if (!is_generic && type_of_values_->is_root_of_generic_subtree()) THROW( ExcGenericWithoutInstance() << EI_Object(type_of_values_->type_name()) );

	return (finished = type_of_values_->finish(is_generic) );
}



string Array::type_name() const {
    return "array_of_" + data_->type_of_values_->type_name();
}



bool Array::operator==(const TypeBase &other) const    {
    return  typeid(*this) == typeid(other) &&
              (*data_->type_of_values_ == static_cast<const Array *>(&other)->get_sub_type() );
}



TypeBase::MakeInstanceReturnType Array::make_instance(std::vector<ParameterPair> vec) const {
	// Create copy of array, we can't set type from parameter vector directly (it's TypeBase that is not allowed)
	Array arr = this->deep_copy();
	// Replace parameter stored in type_of_values_
	MakeInstanceReturnType inst = arr.data_->type_of_values_->make_instance(vec);
	arr.data_->type_of_values_ = inst.first;
	ParameterMap parameter_map = inst.second;
	// Copy attributes
	arr.attributes_ = boost::make_shared<attribute_map>(*attributes_);
	// Set parameters as attribute
	json_string val = this->print_parameter_map_to_json(parameter_map);
	ASSERT( this->validate_json(val), "Invalid JSON format of attribute 'parameters'.\n" );
	(*arr.attributes_)[FlowAttributes::parameters()] = val;
	arr.parameter_map_ = parameter_map;
	(*arr.attributes_)[FlowAttributes::generic_type()] = TypeBase::hash_str();
	arr.generic_type_hash_ = this->content_hash();

	return std::make_pair( boost::make_shared<Array>(arr), parameter_map );
}


Array Array::deep_copy() const {
	Array arr = Array(Integer()); // Type integer will be overwritten
	arr.data_ = boost::make_shared<Array::ArrayData>(*this->data_);
	arr.data_->finished = false;
	return arr;
}


Array::Array(boost::shared_ptr<TypeBase> type, unsigned int min_size, unsigned int max_size)
: data_(boost::make_shared<ArrayData>(min_size, max_size))
{
    ASSERT( min_size <= max_size, "Wrong limits for size of Input::Type::Array, min: %d, max: %d\n", min_size, max_size);
    ASSERT( type->is_closed(), "Sub-type '%s' of Input::Type::Array must be closed!", type->type_name().c_str());

	data_->type_of_values_ = type;
}


/**********************************************************************************
 * implementation and explicit instantiation of Array constructor template
 */

template <class ValueType>
Array::Array(const ValueType &type, unsigned int min_size, unsigned int max_size)
: Array(boost::static_pointer_cast<TypeBase>( boost::make_shared<ValueType>(type) ), min_size, max_size)
{
    // ASSERT MESSAGE: The type of declared keys has to be a class derived from TypeBase.
    BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, ValueType >::value) );
}

// explicit instantiation

#define ARRAY_CONSTRUCT(TYPE) \
template Array::Array(const TYPE &type, unsigned int min_size, unsigned int max_size)

ARRAY_CONSTRUCT(String);
ARRAY_CONSTRUCT(Integer);
ARRAY_CONSTRUCT(Double);
ARRAY_CONSTRUCT(Bool);
ARRAY_CONSTRUCT(FileName);
ARRAY_CONSTRUCT(Selection);
ARRAY_CONSTRUCT(Array);
ARRAY_CONSTRUCT(Record);
ARRAY_CONSTRUCT(Abstract);
ARRAY_CONSTRUCT(Parameter);
ARRAY_CONSTRUCT(Instance);


/**********************************************************************************
 * implementation of Type::Scalar ... and descendants.
 */

/**********************************************************************************
 * implementation of Type::Bool
 */


TypeBase::TypeHash Bool::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, type_name());
    return seed;
}


string Bool::type_name() const {
    return "Bool";
}


TypeBase::MakeInstanceReturnType Bool::make_instance(std::vector<ParameterPair> vec) const {
	return std::make_pair( boost::make_shared<Bool>(*this), ParameterMap() );
}

/**********************************************************************************
 * implementation of Type::Integer
 */

TypeBase::TypeHash Integer::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, lower_bound_);
    boost::hash_combine(seed, upper_bound_);
    return seed;
}



bool Integer::match(std::int64_t value) const {
    return ( value >=lower_bound_ && value <= upper_bound_);
}



string Integer::type_name() const {
    return "Integer";
}


TypeBase::MakeInstanceReturnType Integer::make_instance(std::vector<ParameterPair> vec) const {
	return std::make_pair( boost::make_shared<Integer>(*this), ParameterMap() );
}


/**********************************************************************************
 * implementation of Type::Double
 */


TypeBase::TypeHash Double::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, lower_bound_);
    boost::hash_combine(seed, upper_bound_);
    return seed;
}



bool Double::match(double value) const {
    return ( value >=lower_bound_ && value <= upper_bound_);
}



string Double::type_name() const {
    return "Double";
}


TypeBase::MakeInstanceReturnType Double::make_instance(std::vector<ParameterPair> vec) const {
	return std::make_pair( boost::make_shared<Double>(*this), ParameterMap() );
}


/**********************************************************************************
 * implementation of Type::FileName
 */

TypeBase::TypeHash FileName::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, type_);
    return seed;
}





string FileName::type_name() const {
    switch (type_) {
    case ::FilePath::input_file:
        return "FileName_input";
    case ::FilePath::output_file:
        return "FileName_output";
    default:
        return "FileName";
    }
}



bool FileName::match(const string &str) const {
    return (type_ == ::FilePath::input_file) || (str[0] != DIR_DELIMITER); // output files can not be absolute
}



TypeBase::MakeInstanceReturnType FileName::make_instance(std::vector<ParameterPair> vec) const {
	return std::make_pair( boost::make_shared<FileName>(*this), ParameterMap() );
}


/**********************************************************************************
 * implementation of Type::String
 */


TypeBase::TypeHash String::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, type_name());
    return seed;
}



string String::type_name() const {
    return "String";
}




bool String::match(const string &str) const {
    return true;
}



TypeBase::MakeInstanceReturnType String::make_instance(std::vector<ParameterPair> vec) const {
	return std::make_pair( boost::make_shared<String>(*this), ParameterMap() );
}



} // closing namespace Type
} // closing namespace Input



