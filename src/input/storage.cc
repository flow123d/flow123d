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
 * @file    storage.cc
 * @brief   
 */

#include "system/global_defs.h"
#include "system/system.hh"
#include "input/storage.hh"

namespace Input {
using namespace std;

/**********************************************************************
 * Implements StorageBase
 */

std::int64_t StorageBase::get_int() const {
    THROW( ExcStorageTypeMismatch() << EI_RequestedType("int") << EI_StoredType( typeid(*this).name()) );
    return 0;
}



double StorageBase::get_double() const {
    THROW( ExcStorageTypeMismatch() << EI_RequestedType("double") << EI_StoredType( typeid(*this).name()) );
    return 0;
}



bool StorageBase::get_bool() const {
    THROW( ExcStorageTypeMismatch() << EI_RequestedType("bool") << EI_StoredType( typeid(*this).name()) );
    return false;
}



const std::string & StorageBase::get_string() const {
    THROW( ExcStorageTypeMismatch() << EI_RequestedType("string") << EI_StoredType( typeid(*this).name()) );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-local-addr"
    return 0;   // Return reference to temporary, but we throw anyway.
#pragma GCC diagnostic pop
}



StorageBase * StorageBase::get_item(const unsigned int index) const {
    THROW( ExcStorageTypeMismatch() << EI_RequestedType("array") << EI_StoredType( typeid(*this).name()) );
    return 0;
}



bool StorageBase::is_null() const {
    return false;
}



unsigned int StorageBase::get_array_size() const {
    THROW( ExcStorageTypeMismatch() << EI_RequestedType("array") << EI_StoredType( typeid(*this).name()) );
    return 0;
}

StorageBase::~StorageBase()
{}

/*****************************************************************
 * Implementation of StorageArray
 */

StorageArray::StorageArray(unsigned int size)
: array_(size)
{
    for( vector<StorageBase *>::iterator it = array_.begin(); it != array_.end(); ++it)
        *it=NULL;
}

 StorageBase * StorageArray::deep_copy() const {
    StorageArray *copy = new StorageArray(this->get_array_size());

    for(unsigned int i=0; i< array_.size(); i++)
        if (array_[i] != NULL) copy->new_item(i, array_[i]->deep_copy() );

    return copy;
}

 void StorageArray::new_item(unsigned int index, StorageBase* item) {
	 FEAL_DEBUG_ASSERT(index < array_.size())(index)(array_.size()).error("Index is out of array.");
	 FEAL_ASSERT( array_[index] == NULL )(index).error("Can not replace non NULL pointer.");
     array_[index] = item;
 }


 void StorageArray::set_item(unsigned int index, StorageBase* item) {
	FEAL_DEBUG_ASSERT(index < array_.size())(index)(array_.size()).error("Index is out of array.");
	FEAL_ASSERT( array_[index] == NULL || typeid(*array_[index]) == typeid(StorageNull) )(index)
			.error("Can not replace non NULL pointer.");

    if ( typeid(*array_[index]) == typeid(StorageNull) ) delete array_[index];
    array_[index] = item;
 }



StorageBase * StorageArray::get_item(const unsigned int index) const {
	FEAL_DEBUG_ASSERT( index < array_.size() )(index)(array_.size()).error("Index is out of array.");
    FEAL_DEBUG_ASSERT(array_[index] != NULL)(index).error();
    return array_[index];
}



unsigned int StorageArray::get_array_size() const {
    return array_.size();
}



bool StorageArray::is_null() const {
    return false;
}


void StorageArray::print(ostream &stream, int pad)  const {
    stream << setw(pad) << "" << "array(" << this->get_array_size() << ")" << std::endl;
    for(unsigned int i=0;i<get_array_size();++i) get_item(i)->print(stream, pad+2);
}



StorageArray::~StorageArray() {
    for( vector<StorageBase *>::iterator it = array_.begin(); it != array_.end(); ++it)
        if (*it != NULL) delete (*it);
}

/**********************************************
 * Implementation of StorageBool
 */

StorageBool::StorageBool(bool value)
: value_(value)
{}



bool StorageBool::get_bool() const {
    return value_;
}



bool StorageBool::is_null() const {
    return false;
}



StorageBase * StorageBool::deep_copy() const {
    return new StorageBool(value_);
}

void StorageBool::print(ostream &stream, int pad) const {
    stream << setw(pad) << "" <<  "bool(" << value_ << ")"<<std::endl;
}

StorageBool::~StorageBool()
{}



/**********************************************
 * Implementation of StorageInt
 */

StorageInt::StorageInt(int value)
: value_(value)
{}



std::int64_t StorageInt::get_int() const {
    return value_;
}



bool StorageInt::is_null() const {
    return false;
}

 StorageBase * StorageInt::deep_copy() const {
    return new StorageInt(value_);
}


void StorageInt::print(ostream &stream, int pad) const {
   stream << setw(pad) << "" <<  "int(" << value_ << ")"<<std::endl;
}


StorageInt::~StorageInt()
{}



/**********************************************
 * Implementation of StorageDouble
 */

StorageDouble::StorageDouble(double value)
: value_(value)
{}



double StorageDouble::get_double() const {
    return value_;
}



bool StorageDouble::is_null() const {
    return false;
}

 StorageBase * StorageDouble::deep_copy() const {
    return new StorageDouble(value_);
}


void StorageDouble::print(ostream &stream, int pad) const {
    stream << setw(pad) << "" <<  "double(" << value_ << ")"<<std::endl;
}


StorageDouble::~StorageDouble()
{}



/**********************************************
 * Implementation of StorageString
 */

StorageString::StorageString(const string& value)
: value_(value)
{}



const string& StorageString::get_string() const {
    return value_;
}



bool StorageString::is_null() const {
    return false;
}


StorageBase * StorageString::deep_copy() const {
    return new StorageString(value_);
}


void StorageString::print(ostream &stream, int pad) const {
    stream << setw(pad) << "" << "string(" << value_ << ")"<<std::endl;
}


StorageString::~StorageString()
{}


/**********************************************
 * Implementation of StorageNull
 */


bool StorageNull::is_null() const {
    return true;
}



 StorageBase * StorageNull::deep_copy() const {
    return new StorageNull();
}



void StorageNull::print(ostream &stream, int pad) const{
    stream << setw(pad) << "" << "null()"<<std::endl;
}


StorageNull::~StorageNull()
{}


} // namespace Input
