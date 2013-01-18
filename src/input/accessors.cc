/*
 * accessors.cc
 *
 *  Created on: Apr 26, 2012
 *      Author: jb
 */



#include "input/accessors.hh"


namespace Input {
/*****************************************************************************
 * Implementation of the class Input::Record
 */


Record::Record()
: record_type_(), storage_( NULL )
{}



Record::Record(const Record &rec)
: record_type_(rec.record_type_), storage_(rec.storage_)
{}



Record::Record(const StorageBase *store, const Type::Record type)
: record_type_(type), storage_(store)
{
    if (store->is_null())
        THROW( ExcAccessorForNullStorage() << EI_AccessorName("Record") );
}




/*****************************************************************************
 * Implementation of the class Input::AbstractRecord
 */

AbstractRecord::AbstractRecord()
: record_type_(), storage_( NULL )
{}



AbstractRecord::AbstractRecord(const AbstractRecord &rec)
: record_type_(rec.record_type_), storage_(rec.storage_)
{}



AbstractRecord::AbstractRecord(const StorageBase *store, const Type::AbstractRecord type)
: record_type_(type), storage_(store)
{
    if (store->is_null())
        THROW( ExcAccessorForNullStorage() << EI_AccessorName("AbstractRecord") );
}



AbstractRecord::operator Record()
{ return Record(storage_,type()); }



Input::Type::Record AbstractRecord::type()
{
    unsigned int type_id = storage_->get_item(0)->get_int();
    return record_type_.get_descendant(type_id);
}



/*****************************************************************************
 * Implementation of the class Input::Array
 */


Array::Array()
: array_type_(Type::Bool()), storage_( NULL )
{}


Array::Array(const Array &ar)
: array_type_(ar.array_type_), storage_(ar.storage_)
{}


Array::Array(const StorageBase *store, const Type::Array type)
: array_type_(type), storage_(store)
{
    if (store->is_null())
        THROW( ExcAccessorForNullStorage() << EI_AccessorName("Array") );
}

/*****************************************************************************
 * Explicit instantiation of accessor's templates
 *
 * .. TODO
 */



} // closing namespace Input
