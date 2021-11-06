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
 * @file    accessors.cc
 * @brief   
 */

#include <memory>
#include <boost/lexical_cast.hpp>
#include "input/accessors.hh"


namespace Input {


/*************************************************************************************************************************************
 * Implementation of InputException
 */


std::ostringstream &Exception::form_message(std::ostringstream &converter) const {

    converter << "--------------------------------------------------------" << std::endl;
    converter << "User Error: ";
    print_info(converter);
#ifdef FLOW123D_DEBUG_MESSAGES
    converter << "\n** Diagnosting info **\n" ;
    converter << boost::diagnostic_information_what( *this );
    print_stacktrace(converter);
#endif
    converter << std::endl << "--------------------------------------------------------" << std::endl;

    return converter;
}





/*****************************************************************************
 * Implementation of the class Input::Address
 */

Address::Address()
: data_(std::make_shared<AddressData>())
{
   data_->root_type_ = nullptr;
   data_->root_storage_ = &Array::empty_storage_;
   data_->descendant_order_ = 0;
   data_->actual_storage_ = &Array::empty_storage_;
}


Address::Address(const StorageBase * storage_root, const Type::TypeBase *type_root)
: data_( std::make_shared<AddressData>() )
{
    if (! storage_root)
        THROW( ExcAddressNullPointer() << EI_AccessorName("storage_root") );
    if (! type_root )
        THROW( ExcAddressNullPointer() << EI_AccessorName("type_root") );

    data_->root_type_ = type_root;
    data_->root_storage_ = storage_root;
    data_->descendant_order_ = 0;
    data_->actual_storage_ = storage_root;
}



Address::Address(const Address& other)
: data_(other.data_)
{}


Address::AddressData::~AddressData() {
	if (	!parent_
			&& root_storage_ == actual_storage_
			&& root_type_ ) {
		delete root_storage_;
	}
}


std::shared_ptr<Address> Address::down(unsigned int idx) const {

	auto addr = std::make_shared<Address>(this->data_->root_storage_, this->data_->root_type_);
	addr->data_->parent_ = this->data_;
	addr->data_->descendant_order_ = idx;
	addr->data_->actual_storage_ = data_->actual_storage_->get_item(idx);

	return addr;
}


std::string Address::make_full_address() const {
	std::vector<unsigned int> path;
	std::shared_ptr<AddressData> address_data = data_;
	while (address_data->parent_ != NULL) {
		path.push_back(address_data->descendant_order_);
		address_data = address_data->parent_;
	}

	// for empty path is returned address of root node
	if (path.size() == 0) {
		return "/";
	}

    const StorageBase * storage = address_data->root_storage_;
    const Type::TypeBase * input_type = address_data->root_type_;
	std::string address = "";
	int i = path.size()-1;

    while (i >= 0) {

    	// dispatch types
        if (typeid(*input_type) == typeid(Type::Record)) {
        	storage = storage->get_item(path[i]);
        	const Type::Record * rec = static_cast<const Type::Record *>(input_type);
        	Type::Record::KeyIter it = rec->begin() + path[i];
        	address = address + "/" + it->key_;
        	input_type = it->type_.get();
        	i--;
        } else
		if (typeid(*input_type) == typeid(Type::Abstract)) {
			const Type::Abstract * a_rec = static_cast<const Type::Abstract *>(input_type);
			const StorageString * storage_type = static_cast<const StorageString *>(storage->get_item(0));
			input_type = & a_rec->get_descendant(storage_type->get_string());
		} else
		if (typeid(*input_type) == typeid(Type::Array)) {
	    	storage = storage->get_item(path[i]);
			const Type::Array * arr = static_cast<const Type::Array *>(input_type);
			address = address + "/" + boost::lexical_cast<std::string>(path[i]);
			input_type = & arr->get_sub_type();
			i--;
		}
    }

    return address;
}



/*****************************************************************************
 * Implementation of the class Input::Record
 */


Record::Record()
: address_( Address() ), record_type_()
{}



Record::Record(const Record &rec)
: address_(rec.address_), record_type_(rec.record_type_)
{}



Record::Record(const Address &address, const Type::Record type)
: address_(address), record_type_(type)
{
    if (address.storage_head()->is_null())
        THROW( ExcAccessorForNullStorage() << EI_AccessorName("Record") );
}


Input::EI_Address Record::ei_address() const
{
	return EI_Address(address_string());
}


string Record::address_string() const
{
	return address_.make_full_address();
}

string Record::input_type_name()
{
	return record_type_.type_name();
}


Type::Record::KeyIter Record::get_type_key_iterator(const string &key) const {
	return record_type_.key_iterator(key);
}



/*****************************************************************************
 * Implementation of the class Input::Tuple
 */


Tuple::Tuple()
: tuple_type_()
{
	this->address_ = Address();
}



Tuple::Tuple(const Tuple &tpl)
: Record(tpl), tuple_type_(tpl.tuple_type_)
{
	this->address_ = tpl.address_;
}



Tuple::Tuple(const Address &address, const Type::Tuple type)
: tuple_type_(type)
{
	this->address_ = address;
    if (address.storage_head()->is_null())
        THROW( ExcAccessorForNullStorage() << EI_AccessorName("Tuple") );
}


string Tuple::input_type_name()
{
	return tuple_type_.type_name();
}


Type::Record::KeyIter Tuple::get_type_key_iterator(const string &key) const {
	return tuple_type_.key_iterator(key);
}



/*****************************************************************************
 * Implementation of the class Input::AbstractRecord
 */

AbstractRecord::AbstractRecord()
: abstract_type_(), address_( Address() )
{}



AbstractRecord::AbstractRecord(const AbstractRecord &rec)
: abstract_type_(rec.abstract_type_), address_(rec.address_)
{}



AbstractRecord::AbstractRecord(const Address &address, const Type::Abstract type)
: abstract_type_(type), address_(address)
{
	if (address.storage_head()->is_null())
        THROW( ExcAccessorForNullStorage() << EI_AccessorName("AbstractRecord") );
}



AbstractRecord::operator Record() const
{ return Record(address_,type()); }



Input::Type::Record AbstractRecord::type() const
{
	string type_name = address_.storage_head()->get_item(0)->get_string();
    return abstract_type_.get_descendant(type_name);
}


Input::EI_Address AbstractRecord::ei_address() const
{
	return EI_Address(address_string());
}

string AbstractRecord::address_string() const
{
	return address_.make_full_address();
}


/*****************************************************************************
 * Implementation of the class Input::Array
 */


Array::Array()
: array_type_(Type::Bool()), address_( Address() )
{}


Array::Array(const Array &ar)
: array_type_(ar.array_type_), address_(ar.address_)
{}


Array::Array(const Address &address, const Type::Array type)
: array_type_(type), address_(address)
{
    if (address.storage_head()->is_null())
        THROW( ExcAccessorForNullStorage() << EI_AccessorName("Array") );
}


Input::EI_Address Array::ei_address() const
{
	return EI_Address(address_string());
}



string Array::address_string() const
{
	return address_.make_full_address();
}


StorageArray Array::empty_storage_ = StorageArray(0);



/*****************************************************************************
 * Implementation of the class Input::IteratorBase
 */


IteratorBase::IteratorBase(const Address &address, const unsigned int index)
: address_(address), index_(index)
{}


const Address &IteratorBase::get_address() const {
	return address_;
}



/*****************************************************************************
 * Explicit instantiation of accessor's templates
 *
 * .. TODO
 */



} // closing namespace Input
