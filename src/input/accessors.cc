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
   data_->root_storage_ = &Array::empty_storage_;
   data_->parent_ = NULL;
   data_->descendant_order_ = 0;
   data_->actual_storage_ = &Array::empty_storage_;
}


Address::Address(const StorageBase * storage_root, const Type::TypeBase *type_root)
: data_( boost::make_shared<AddressData>() )
{
    if (storage_root == NULL)
        THROW( ExcAddressNullPointer() << EI_AccessorName("storage_root") );
    if (!type_root || type_root == NULL)
        THROW( ExcAddressNullPointer() << EI_AccessorName("type_root") );

    data_->root_type_ = type_root;
    data_->root_storage_ = storage_root;
    data_->parent_ = NULL;
    data_->descendant_order_ = 0;
    data_->actual_storage_ = storage_root;
}



Address::Address(const Address& other)
: data_(other.data_)
{}


const Address * Address::down(unsigned int idx) const {

	Address *addr = new Address(this->data_->root_storage_, this->data_->root_type_);
	addr->data_->parent_ = this->data_.get();
	addr->data_->descendant_order_ = idx;
	addr->data_->actual_storage_ = data_->actual_storage_->get_item(idx);

	return addr;
}


std::string Address::make_full_address() const {
	std::vector<unsigned int> path;
	AddressData * address_data = data_.get();
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
		if (typeid(*input_type) == typeid(Type::AbstractRecord)) {
			const Type::AbstractRecord * a_rec = static_cast<const Type::AbstractRecord *>(input_type);
			const StorageInt * storage_type = static_cast<const StorageInt *>(storage->get_item(0));
			input_type = & a_rec->get_descendant(storage_type->get_int());
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


const Input::Address & Record::get_address() const
{
	return address_;
}


void Record::set_address(const Address &address)
{
	address_ = address;
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


const Input::Address & AbstractRecord::get_address() const
{
	return address_;
}


void AbstractRecord::set_address(const Address &address)
{
	address_ = address;
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


const Input::Address & Array::get_address() const
{
	return address_;
}


void Array::set_address(const Address &address)
{
	address_ = address;
}


StorageArray Array::empty_storage_ = StorageArray(0);


/*****************************************************************************
 * Explicit instantiation of accessor's templates
 *
 * .. TODO
 */



} // closing namespace Input
