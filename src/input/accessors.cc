/*
 * accessors.cc
 *
 *  Created on: Apr 26, 2012
 *      Author: jb
 */


#include <memory>
#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>
#include "input/accessors.hh"
#include "input/storage_transpose.hh"


namespace Input {


/*************************************************************************************************************************************
 * Implementation of InputException
 */


const char * Exception::what() const throw () {
    // have preallocated some space for error message we want to return
    // Is there any difference, if we move this into ExceptionBase ??
    static std::string message(1024,' ');


    // Be sure that this function do not throw.
    try {
        std::ostringstream converter;

        converter << std::endl << std::endl;
        converter << "--------------------------------------------------------" << std::endl;
        converter << "User Error: ";
        print_info(converter);
#ifdef FLOW123D_DEBUG_MESSAGES
        converter << "\n** Diagnosting info **\n" ;
        converter << boost::diagnostic_information_what( *this );
        print_stacktrace(converter);
#endif
        converter << "--------------------------------------------------------" << std::endl;

        message = converter.str();
        return message.c_str();

    } catch (std::exception &exc) {
        std::cerr << "*** Exception encountered in exception handling routines ***" << std::endl << "*** Message is " << std::endl
                << exc.what() << std::endl << "*** Aborting! ***" << std::endl;

        std::abort();
    } catch (...) {
        std::cerr << "*** Exception encountered in exception handling routines ***" << std::endl << "*** Aborting! ***"
                << std::endl;

        std::abort();
    }
    return 0;
}





/*****************************************************************************
 * Implementation of the class Input::Address
 */

Address::Address()
: data_(boost::make_shared<AddressData>())
{
   data_->root_type_ = nullptr;
   data_->root_storage_ = &Array::empty_storage_;
   data_->parent_ = nullptr;
   data_->descendant_order_ = 0;
   data_->actual_storage_ = &Array::empty_storage_;
}


Address::Address(const StorageBase * storage_root, const Type::TypeBase *type_root)
: data_( boost::make_shared<AddressData>() )
{
    if (! storage_root)
        THROW( ExcAddressNullPointer() << EI_AccessorName("storage_root") );
    if (! type_root )
        THROW( ExcAddressNullPointer() << EI_AccessorName("type_root") );

    data_->root_type_ = type_root;
    data_->root_storage_ = storage_root;
    data_->parent_ = nullptr;
    data_->descendant_order_ = 0;
    data_->actual_storage_ = storage_root;
}



Address::Address(const Address& other)
: data_(other.data_)
{}


std::shared_ptr<Address> Address::down(unsigned int idx) const {

	auto addr = std::make_shared<Address>(this->data_->root_storage_, this->data_->root_type_);
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


Input::EI_Address Record::ei_address() const
{
	return EI_Address(address_string());
}


string Record::address_string() const
{
	return address_.make_full_address();
}

string Record::record_type_name()
{
	return record_type_.type_name();
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
	string type_name = address_.storage_head()->get_item(0)->get_string();
    return record_type_.get_descendant(type_name);
}


Input::EI_Address AbstractRecord::ei_address() const
{
	return EI_Address(address_string());
}

string AbstractRecord::address_string() const
{
	return address_.make_full_address();
}

void AbstractRecord::transpose_to(Input::Record &target_rec, string target_key, unsigned int vec_size) {
	Input::Iterator<Array> it = target_rec.find<Array>(target_key);
	if (it) { // target_key is set by user
		return;
	}

	Type::Record::KeyIter key_it = target_rec.record_type_.key_iterator(target_key);
	const Type::Array *in_arr = static_cast<const Type::Array *>(key_it->type_.get());
	const Type::TypeBase *target_type = &(in_arr->get_sub_type());

	StorageTranspose trans(target_type, &(this->record_type_), this->address_.storage_head(), vec_size);
    StorageArray* result_storage = new StorageArray(vec_size);
    for(unsigned int i=0; i<vec_size; i++) {
    	result_storage->new_item(i, trans.get_item(i));
    }
    StorageArray* storage_arr =
            const_cast<StorageArray *>(
            static_cast<const StorageArray *>(target_rec.address_.storage_head()));
    storage_arr->set_item(target_rec.record_type_.key_index(target_key), result_storage);
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
 * Explicit instantiation of accessor's templates
 *
 * .. TODO
 */



} // closing namespace Input
