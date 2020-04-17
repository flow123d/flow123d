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
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include <boost/pointer_cast.hpp>


#include "input_type.hh"
#include "type_output.hh"
#include "type_repository.hh"
#include "attribute_lib.hh"
#include "json_spirit/json_spirit.h"

#include <stdint.h>                                          // for int64_t
#include <boost/algorithm/string.hpp>
#include <boost/exception/detail/error_info_impl.hpp>        // for error_info
#include <boost/exception/info.hpp>                          // for operator<<
#include <boost/functional/hash/hash.hpp>                    // for hash_com...
#include <boost/static_assert.hpp>                           // for BOOST_ST...
#include <boost/type_traits/is_base_of.hpp>                  // for is_base_of
#include <memory>                                            // for shared_ptr
#include <ostream>                                           // for operator<<
#include <typeinfo>                                          // for type_info
#include <utility>                                           // for make_pair
#include "input/json_spirit/json_spirit_error_position.h"    // for Error_po...
#include "input/json_spirit/json_spirit_reader.h"            // for read_or_...
#include "input/json_spirit/json_spirit_value.h"             // for mValue
#include "input/type_base.hh"                                // for TypeBase
#include "input/type_generic.hh"                             // for ExcGener...
#include "system/asserts.hh"                                 // for Assert
#include "system/exceptions.hh"                              // for ExcGener...
#include "system/file_path.hh"                               // for FilePath
namespace Input { namespace Type { class Abstract; } }
namespace Input { namespace Type { class Record; } }
namespace Input { namespace Type { class Selection; } }
namespace Input { namespace Type { class Tuple; } }



namespace Input {
namespace Type {

using namespace std;



/*******************************************************************
 * implementation of TypeBase
 */



TypeBase::TypeBase()
: attributes_( std::make_shared<attribute_map>() ), root_of_generic_subtree_(false),
  generic_type_hash_(TypeBase::none_hash) {}



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



void TypeBase::delete_unfinished_types() {
	// mark unfinished types as deleted
	Input::TypeRepository<Instance>::get_instance().finish(FinishStatus::delete_);
	Input::TypeRepository<Abstract>::get_instance().finish(FinishStatus::delete_);
	Input::TypeRepository<Record>::get_instance().finish(FinishStatus::delete_);
	Input::TypeRepository<Tuple>::get_instance().finish(FinishStatus::delete_);
	Input::TypeRepository<Selection>::get_instance().finish(FinishStatus::delete_);

	// check and remove deleted types
	Input::TypeRepository<Instance>::get_instance().reset_deleted_types();
	Input::TypeRepository<Abstract>::get_instance().reset_deleted_types();
	Input::TypeRepository<Record>::get_instance().reset_deleted_types();
	Input::TypeRepository<Tuple>::get_instance().reset_deleted_types();
	Input::TypeRepository<Selection>::get_instance().reset_deleted_types();
}



 std::string TypeBase::hash_str(TypeHash hash) {
    stringstream ss;
    ss << "\"" << std::hex << hash << "\"";
    return ss.str();
}




void TypeBase::add_attribute_(std::string name, json_string val) {
	ASSERT(validate_json(val))(name)(val).error("Invalid JSON format of attribute");
	(*attributes_)[name] = val;
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
	ss << "{";
	for (ParameterMap::iterator it=parameter_map.begin(); it!=parameter_map.end(); it++) {
		if (it != parameter_map.begin()) ss << "," << endl;
		ss << "\"" << (*it).first << "\" : " << TypeBase::hash_str( (*it).second );
	}
	ss << "}";
	return ss.str();
}

TypeBase::json_string TypeBase::print_parameter_map_keys_to_json(ParameterMap parameter_map) const {
    stringstream ss;
    ss << "[ ";
    for (ParameterMap::iterator it=parameter_map.begin(); it!=parameter_map.end(); it++) {
        if (it != parameter_map.begin()) ss << ", ";
        ss << "\"" << it->first << "\"";
    }
    ss << " ]";
    return ss.str();
}

void TypeBase::set_generic_attributes(ParameterMap parameter_map) {
    // check if the type is really generic (it may be non-generic even if part of a generic subtree)
    if (parameter_map.size() > 0)
        add_attribute_(Attribute::generic_parameters(), print_parameter_map_keys_to_json(parameter_map));
    if (is_root_of_generic_subtree())
        add_attribute_(Attribute::root_of_generic_subtree(), "true");
}


void TypeBase::copy_attributes(attribute_map other_attributes) {
    attributes_->clear();
    for(auto &item : other_attributes) {
        if (item.first[0] != '_') // not internal attribute
            attributes_->insert(item);
    }
}


FinishStatus TypeBase::finish_status() const {
	return FinishStatus::regular_;
}


bool TypeBase::is_finished() const {
	return true;
}


bool TypeBase::is_closed() const {
	return true;
}


string TypeBase::type_name() const {
	return "TypeBase";
}


string TypeBase::class_name() const {
	return "TypeBase";
}


FinishStatus TypeBase::finish(FinishStatus finish_type) {
	ASSERT((finish_type != FinishStatus::none_) && (finish_type != FinishStatus::in_perform_)).error();
	return finish_type;
}


bool TypeBase::operator==(const TypeBase &other) const {
	return typeid(*this) == typeid(other);
}

bool TypeBase::operator!=(const TypeBase & other) const {
	return ! (*this == other);
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


FinishStatus Array::finish(FinishStatus finish_type) {
	return data_->finish(finish_type);
}



Array::ArrayData::ArrayData(unsigned int min_size, unsigned int max_size)
: lower_bound_(min_size), upper_bound_(max_size), finish_status(FinishStatus::none_)
{}


FinishStatus Array::ArrayData::finish(FinishStatus finish_type)
{
	ASSERT(finish_type != FinishStatus::none_).error();
	ASSERT(finish_type != FinishStatus::in_perform_).error();
	ASSERT(finish_status != FinishStatus::in_perform_).error("Recursion in the IST element: array_of_" + type_of_values_->type_name());

	if (finish_status != FinishStatus::none_) return finish_status;


	finish_status = FinishStatus::in_perform_;

	if (typeid( *(type_of_values_.get()) ) == typeid(Instance)) {
		type_of_values_->finish(FinishStatus::generic_); // finish Instance object
		type_of_values_ = type_of_values_->make_instance().first;
	}
	if ((finish_type != FinishStatus::generic_) && type_of_values_->is_root_of_generic_subtree())
		THROW( ExcGenericWithoutInstance() << EI_Object(type_of_values_->type_name()) );

	type_of_values_->finish(finish_type);
	ASSERT(type_of_values_->is_finished()).error();
	if (finish_type == FinishStatus::delete_) type_of_values_.reset();
	finish_status = finish_type;
	return (finish_status);
}



string Array::type_name() const {
    return "array_of_" + data_->type_of_values_->type_name();
}



string Array::class_name() const {
	return "Array";
}



bool Array::operator==(const TypeBase &other) const    {
    return  typeid(*this) == typeid(other) &&
              (*data_->type_of_values_ == static_cast<const Array *>(&other)->get_sub_type() );
}



TypeBase::MakeInstanceReturnType Array::make_instance(std::vector<ParameterPair> vec)  {
	// Create copy of array, we can't set type from parameter vector directly (it's TypeBase that is not allowed)
	Array arr = this->deep_copy();
	// Replace parameter stored in type_of_values_
	MakeInstanceReturnType inst = arr.data_->type_of_values_->make_instance(vec);
	arr.data_->type_of_values_ = inst.first;
	ParameterMap parameter_map = inst.second;
	// Copy attributes
	arr.copy_attributes(*attributes_);

	// Set parameters as attribute
	json_string val = this->print_parameter_map_to_json(parameter_map);
	ASSERT(this->validate_json(val))(val).error("Invalid JSON format of attribute 'parameters'.");
	arr.parameter_map_ = parameter_map;
	arr.generic_type_hash_ = this->content_hash();

	return std::make_pair( std::make_shared<Array>(arr), parameter_map );
}


Array Array::deep_copy() const {
	Array arr = Array(Integer()); // Type integer will be overwritten
	arr.data_ = std::make_shared<Array::ArrayData>(*this->data_);
	arr.data_->finish_status = FinishStatus::none_;
	return arr;
}


Array::Array(std::shared_ptr<TypeBase> type, unsigned int min_size, unsigned int max_size)
: data_(std::make_shared<ArrayData>(min_size, max_size))
{
	ASSERT_LE(min_size, max_size).error("Wrong limits for size of Input::Type::Array");
	ASSERT(type->is_closed()).error();

	data_->type_of_values_ = type;
}


/// Override @p Type::TypeBase::finish_status.
FinishStatus Array::finish_status() const {
    return data_->finish_status; }


bool Array::is_finished() const {
	return (data_->finish_status != FinishStatus::none_) && (data_->finish_status != FinishStatus::in_perform_);
}


/**********************************************************************************
 * implementation and explicit instantiation of Array constructor template
 */

template <class ValueType>
Array::Array(const ValueType &type, unsigned int min_size, unsigned int max_size)
: Array(std::static_pointer_cast<TypeBase>( std::make_shared<ValueType>(type) ), min_size, max_size)
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
ARRAY_CONSTRUCT(Tuple);
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


string Bool::class_name() const {
	return "Bool";
}


TypeBase::MakeInstanceReturnType Bool::make_instance(FMT_UNUSED std::vector<ParameterPair> vec)  {
	return std::make_pair( std::make_shared<Bool>(*this), ParameterMap() );
}

/**********************************************************************************
 * implementation of Type::Integer
 */

Integer::Integer(int lower_bound, int upper_bound)
: lower_bound_(lower_bound), upper_bound_(upper_bound)
{}



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


string Integer::class_name() const {
	return "Integer";
}


TypeBase::MakeInstanceReturnType Integer::make_instance(FMT_UNUSED std::vector<ParameterPair> vec) {
	return std::make_pair( std::make_shared<Integer>(*this), ParameterMap() );
}


/**********************************************************************************
 * implementation of Type::Double
 */


Double::Double(double lower_bound, double upper_bound)
: lower_bound_(lower_bound), upper_bound_(upper_bound)
{}


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


string Double::class_name() const {
	return "Double";
}


TypeBase::MakeInstanceReturnType Double::make_instance(FMT_UNUSED std::vector<ParameterPair> vec) {
	return std::make_pair( std::make_shared<Double>(*this), ParameterMap() );
}


/**********************************************************************************
 * implementation of Type::FileName
 */

FileName::FileName() {}



FileName::FileName(enum ::FilePath::FileType type)
: type_(type)
{}



FileName FileName::input()
{
	return FileName(::FilePath::input_file);
}



FileName FileName::output()
{
	return FileName(::FilePath::output_file);
}



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


string FileName::class_name() const {
	return "FileName";
}



bool FileName::match(const string &str) const {
    return (type_ == ::FilePath::input_file) || (str[0] != DIR_DELIMITER); // output files can not be absolute
}



TypeBase::MakeInstanceReturnType FileName::make_instance(FMT_UNUSED std::vector<ParameterPair> vec)  {
	return std::make_pair( std::make_shared<FileName>(*this), ParameterMap() );
}



::FilePath::FileType FileName::get_file_type() const {
    return type_;
}




bool FileName::operator==(const TypeBase &other) const
{
	return typeid(*this) == typeid(other) &&
                 (type_== static_cast<const FileName *>(&other)->get_file_type() );
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



string String::class_name() const {
	return "String";
}



bool String::match(FMT_UNUSED const string &str) const {
    return true;
}



TypeBase::MakeInstanceReturnType String::make_instance(FMT_UNUSED std::vector<ParameterPair> vec) {
	return std::make_pair( std::make_shared<String>(*this), ParameterMap() );
}



} // closing namespace Type
} // closing namespace Input



