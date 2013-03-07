/*
 * accessors.cc
 *
 *  Created on: Apr 26, 2012
 *      Author: jb
 */


#include <boost/make_shared.hpp>
#include "input/accessors.hh"


namespace Input {

/*****************************************************************************
 * Implementation of the class Input::Address
 */

Address::Address(const StorageBase * storage_root, const Type::TypeBase *type_root)
: data_( boost::make_shared<AddressData>() )
{
   data_->root_type_ = type_root;
   data_->root_storage_ =storage_root;
   actual_storage_ = storage_root;
   actual_node_=0;
}



Address::Address(const Address& other)
: data_(other.data_),
  actual_node_( other.actual_node_),
  actual_storage_( other.actual_storage_)
{}


void Address::down(unsigned int idx) {
    actual_storage_ = actual_storage_->get_item(idx);
    actual_node_++;
    data_->path_.push_back(idx);
}


std::string Address::make_full_address() {
    std::string address = "/";
    //dodelat

    return address;
}


/*****************************************************************************
 * Implementation of the class Input::Record
 */


Record::Record()
: record_type_(), address_( Array::empty_address_ )
{}



Record::Record(const Record &rec)
: record_type_(rec.record_type_), address_(rec.address_)
{}



Record::Record(const Address &address, const Type::Record type)
: record_type_(type), address_(address)
{
    if (address.storage_head()->is_null())
        THROW( ExcAccessorForNullStorage() << EI_AccessorName("Record") );
}




/*****************************************************************************
 * Implementation of the class Input::AbstractRecord
 */

AbstractRecord::AbstractRecord()
: record_type_(), address_( Array::empty_address_ )
{}



AbstractRecord::AbstractRecord(const AbstractRecord &rec)
: record_type_(rec.record_type_), address_(rec.address_)
{}



AbstractRecord::AbstractRecord(const Address &address, const Type::AbstractRecord type)
: record_type_(type), address_(address)
{
    if (address.storage_head()->is_null())
        THROW( ExcAccessorForNullStorage() << EI_AccessorName("AbstractRecord") );
}



AbstractRecord::operator Record() const
{ return Record(address_,type()); }



Input::Type::Record AbstractRecord::type() const
{
    unsigned int type_id = address_.storage_head()->get_item(0)->get_int();
    return record_type_.get_descendant(type_id);
}



/*****************************************************************************
 * Implementation of the class Input::Array
 */


Array::Array()
: array_type_(Type::Bool()), address_( empty_address_ )
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

StorageArray Array::empty_storage_ = StorageArray(0);

Address Array::empty_address_ = Address( &Array::empty_storage_, NULL );






/*****************************************************************************
 * Explicit instantiation of accessor's templates
 *
 * .. TODO
 */



} // closing namespace Input
