/*
 * storage.cc
 *
 *  Created on: May 5, 2012
 *      Author: jb
 */

#include "storage.hh"

namespace Input {
using namespace std;

/**********************************************************************
 * Implements StorageBase
 */

int StorageBase::get_int() const {
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
    return 0;
}



const StorageBase * StorageBase::get_item(const unsigned int index) const {
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

 StorageBase * StorageArray::deep_copy() {
    StorageArray *copy = new StorageArray(this->get_array_size());

    for(int i=0; i< array_.size(); i++)
        if (array_[i] != NULL) copy->new_item(i, array_[i]->deep_copy() );

    return copy;
}

void StorageArray::new_item(unsigned int index, StorageBase* item) {
    ASSERT( index < array_.size() , "Index %d out of array of size: %d", index, array_.size());
    if (array_[index] == NULL) array_[index] = item;
    else xprintf(PrgErr, "Can not replace non NULL pointer.");
}



const StorageBase * StorageArray::get_item(const unsigned int index) const {
    ASSERT( index < array_.size() , "Index %d out of array of size: %d", index, array_.size());
    return array_[index];
}



unsigned int StorageArray::get_array_size() const {
    return array_.size();
}



bool StorageArray::is_null() const {
    return false;
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



 StorageBase * StorageBool::deep_copy() {
    return new StorageBool(value_);
}


StorageBool::~StorageBool()
{}



/**********************************************
 * Implementation of StorageInt
 */

StorageInt::StorageInt(int value)
: value_(value)
{}



int StorageInt::get_int() const {
    return value_;
}



bool StorageInt::is_null() const {
    return false;
}

 StorageBase * StorageInt::deep_copy() {
    return new StorageInt(value_);
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

 StorageBase * StorageDouble::deep_copy() {
    return new StorageDouble(value_);
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


 StorageBase * StorageString::deep_copy() {
    return new StorageString(value_);
}


StorageString::~StorageString()
{}


/**********************************************
 * Implementation of StorageNull
 */


bool StorageNull::is_null() const {
    return true;
}



 StorageBase * StorageNull::deep_copy() {
    return new StorageNull();
}



StorageNull::~StorageNull()
{}


} // namespace Input
