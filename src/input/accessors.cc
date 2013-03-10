/*
 * accessors.cc
 *
 *  Created on: Apr 26, 2012
 *      Author: jb
 */


#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>
#include "input/accessors.hh"


namespace Input {

/*****************************************************************************
 * Implementation of the class Input::Address
 */

Address::Address()
: data_(boost::make_shared<AddressData>())
{
   data_->root_type_ = NULL;
   data_->root_storage_ = new StorageArray(0);
   actual_storage_ = new StorageArray(0);
   actual_node_ = 0;
}


Address::Address(const StorageBase * storage_root, const Type::TypeBase *type_root)
: data_( boost::make_shared<AddressData>() )
{
    data_->root_type_ = type_root;
    data_->root_storage_ = storage_root;
    actual_storage_ = storage_root;
    actual_node_ = 0;
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
    std::string address = "";
    const StorageBase * storage = data_->root_storage_;
    const Type::TypeBase * input_type = data_->root_type_;

    for (int i = 0; i < actual_node_; i++) {
    	storage = storage->get_item(data_->path_[i]);

    	// dispatch types
        if (typeid(*input_type) == typeid(Type::Record)) {
        	const Type::Record * rec = static_cast<const Type::Record *>(input_type);
        	Type::Record::KeyIter it = rec->begin() + data_->path_[i];
        	address = address + "/" + it->key_;
        	input_type = it->type_.get();
        } else
		if (typeid(*input_type) == typeid(Type::AbstractRecord)) {
			const Type::AbstractRecord * a_rec = static_cast<const Type::AbstractRecord *>(input_type);
			const StorageInt * storage_type = static_cast<const StorageInt *>(storage->get_item(0));
			input_type = & a_rec->get_descendant(storage_type->get_int());
		} else
		if (typeid(*input_type) == typeid(Type::Array)) {
			const Type::Array * arr = static_cast<const Type::Array *>(input_type);
			address = address + "/" + boost::lexical_cast<std::string>(data_->path_[i]);
			input_type = & arr->get_sub_type();
		}
    }

    return address;
}



/*****************************************************************************
 * Implementation of the class Input::Record
 */


Record::Record()
: record_type_(), address_( Address() )
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
: record_type_(), address_( Address() )
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




/*****************************************************************************
 * Explicit instantiation of accessor's templates
 *
 * .. TODO
 */



} // closing namespace Input
