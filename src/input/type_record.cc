/*
 * input_type.cc
 *
 *  Created on: Mar 29, 2012
 *      Author: jb
 */


#include "input_type.hh"
#include "type_repository.hh"

#include <limits>
#include <ios>
#include <map>
#include <vector>
#include <string>
#include <iomanip>

#include "system/system.hh"

#include <boost/typeof/typeof.hpp>
#include <boost/type_traits.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>

namespace Input {
namespace Type {

using namespace std;

/*******************************************************************
 * implementation of Default
 */

Default::Default()
: value_("OPTIONAL"), type_(no_default_optional_type)
{}

Default::Default(const std::string & value)
: value_(value), type_(default_at_declaration)
{
    boost::algorithm::trim(value_);
}



Default::Default(enum DefaultType type, const std::string & value)
: value_(value), type_(type)
{}

TypeBase::TypeHash Default::content_hash() const
{   TypeBase::TypeHash seed = 0;
    boost::hash_combine(seed, "Default");
    boost::hash_combine(seed, type_);
    boost::hash_combine(seed, value_);
    return seed;
}



/**********************************************************************************
 * implementation of Type::Record
 */


Record::Record()
: data_( boost::make_shared<RecordData> ("EmptyRecord","") )
{
	close();
    finish();
}



Record::Record(const Record & other)
: TypeBase( other ), data_(other.data_)
{}



Record::Record(const string & type_name_in, const string & description)
: data_( boost::make_shared<RecordData>(type_name_in, description) )

{}


TypeBase::TypeHash Record::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, "Record");
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, data_->description_);
    boost::hash_combine(seed, data_->auto_conversion_key);
    for( Key &key : data_->keys) {
        boost::hash_combine(seed, key.key_);
        boost::hash_combine(seed, key.description_);
        boost::hash_combine(seed, key.type_->content_hash() );
        boost::hash_combine(seed, key.default_.content_hash() );

    }
    return seed;
}



Record &Record::allow_auto_conversion(const string &from_key) {
    ASSERT(data_->auto_conversion_key_idx == -1, "Can not use key %s for auto conversion, the key is already set.", from_key.c_str());
    data_->auto_conversion_key_idx = 0;
    data_->auto_conversion_key=from_key;

    return *this;
}



void Record::make_derive_from(AbstractRecord &parent) {
    if (data_->derived_) return;

    parent.finish();

    // add TYPE key derived from parent
    KeyHash key_h = key_hash("TYPE");
    std::vector<Key>::iterator it = data_->keys.begin();
    RecordData::key_to_index_const_iter kit = data_->key_to_index.find(key_h);
    if (kit != data_->key_to_index.end()) { // check if TYPE key exists and remove its
    	data_->keys.erase( data_->keys.begin()+kit->second );
    }
    Key type_key = { 0, "TYPE", "Sub-record selection.", parent.child_data_->selection_of_childs, Default( type_name() ), true };
    data_->key_to_index[key_h] = type_key.key_index;
    data_->keys.insert(it, type_key);
    for (unsigned int i=0; i<data_->keys.size(); i++) {
		data_->keys[i].key_index = i;
		data_->key_to_index[key_hash( data_->keys[i].key_)] = i;
    }

    data_->derived_ = true;
}



void Record::make_copy_keys(Record &origin) {

    ASSERT(origin.is_closed(), "Origin record is not closed!\n");

	std::vector<Key>::iterator it = data_->keys.begin();
	int n_inserted = 0;
	for(KeyIter pit=origin.data_->keys.begin(); pit != origin.data_->keys.end(); ++pit) {
		Key tmp_key=*pit;    // make temporary copy of the key
		KeyHash key_h = key_hash(tmp_key.key_);

		tmp_key.derived = true;

		// we have to copy TYPE also since there should be place in storage for it
		// however we change its Default to name of actual Record
		if (tmp_key.key_=="TYPE")
			tmp_key.default_=Default( type_name() );

		// check for duplicate keys, override keys of the parent record by the child record
		RecordData::key_to_index_const_iter kit = data_->key_to_index.find(key_h);
		if (kit != data_->key_to_index.end()) {
			// in actual record exists a key with same name as in parent record
			// use values form the child record
			Key *k = &(data_->keys[kit->second+n_inserted]); // indices in key_to_index are not yet updated

			tmp_key.key_ = k->key_;
			tmp_key.description_ = k->description_;
			tmp_key.type_ = k->type_;
			tmp_key.default_ = k->default_;
			tmp_key.derived = false;
			k->key_ = ""; // mark original key for deletion
		}

		data_->key_to_index[key_h] = tmp_key.key_index;

		it = data_->keys.insert(it, tmp_key)+1;
		n_inserted++;
	}
	// delete duplicate keys and update key indices
	for (unsigned int i=0; i<data_->keys.size(); i++) {
		if (data_->keys[i].key_.compare("") == 0) {
			data_->keys.erase( data_->keys.begin()+i);
			i--;
		} else {
			data_->keys[i].key_index = i;
			data_->key_to_index[key_hash( data_->keys[i].key_)] = i;
		}
	}
}



Record &Record::derive_from(AbstractRecord &parent) {
	ASSERT( ! data_->parent_ptr_ , "Record has been already derived.\n");
	ASSERT( parent.is_closed(), "Parent AbstractRecord '%s' must be closed!\n", parent.type_name().c_str());
    data_->parent_ptr_=boost::make_shared<AbstractRecord>(parent);

	return *this;
}



Record &Record::copy_keys(const Record &other) {
	ASSERT( other.is_closed(), "Record '%s' must be closed!\n", other.type_name().c_str());

   	Record tmp(other);
   	make_copy_keys(tmp);

   	return *this;
}


bool Record::is_finished() const {
    return data_->finished;
}


bool Record::is_closed() const {
	return data_->closed_;
}




bool Record::check_key_default_value(const Default &dflt, const TypeBase &type, const string & k_name) const
{
    if ( dflt.has_value_at_declaration() ) {

        try {
            return type.valid_default( dflt.value() );
        } catch (ExcWrongDefault & e) {
            e << EI_KeyName(k_name);
            throw;
        }
    }

    return false;
}




bool Record::finish()
{
	if (data_->finished) return true;

	ASSERT(data_->closed_, "Finished Record '%s' must be closed!", this->type_name().c_str());

    // finish inheritance if parent is non-null
    if (data_->parent_ptr_) make_derive_from(* (data_->parent_ptr_));

    data_->finished = true;
    for (vector<Key>::iterator it=data_->keys.begin(); it!=data_->keys.end(); it++)
    {
        if (it->key_ != "TYPE") {
            data_->finished = data_->finished && const_cast<TypeBase *>( it->type_.get() )->finish();

            // we check once more even keys that was already checked, otherwise we have to store
            // result of validity check in every key
            check_key_default_value(it->default_, *(it->type_), it->key_);
        }
    }

    // Check default values of autoconvertible records
    if (data_->auto_conversion_key_idx != -1 ) {
        data_->auto_conversion_key_idx=key_index(data_->auto_conversion_key);

        // check that all other obligatory keys have default values
        for(KeyIter it=data_->keys.begin(); it != data_->keys.end(); ++it) {
            if (it->default_.is_obligatory() && (int)(it->key_index) != data_->auto_conversion_key_idx)
                xprintf(PrgErr, "Finishing Record auto convertible from the key '%s', but other obligatory key: '%s' has no default value.\n",
                        data_->auto_conversion_key_iter()->key_.c_str(), it->key_.c_str());
        }
    }

    return (data_->finished);
}



const Record &Record::close() const {
    data_->closed_=true;
    return *( Input::TypeRepository<Record>::getInstance().add_type( *this ).get() );
}



string Record::type_name() const {
    return data_->type_name_;
}


string Record::full_type_name() const {
	if (data_->parent_ptr_) {
		return data_->type_name_ + ":" + data_->parent_ptr_->type_name();
	}
    return data_->type_name_;
}



bool Record::valid_default(const string &str) const
{
    if (data_->auto_conversion_key_idx >=0) {
        unsigned int idx=key_index(data_->auto_conversion_key);
        if ( data_->keys[idx].type_ ) return data_->keys[idx].type_->valid_default(str);
    } else {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(this->type_name()));
    }
    return false;
}



bool Record::operator==(const TypeBase &other) const
{ return  typeid(*this) == typeid(other) &&
                 (type_name() == static_cast<const Record *>(&other)->type_name() );
}


Record::KeyIter Record::auto_conversion_key_iter() const {
    finished_check();
    return data_->auto_conversion_key_iter();
}


Record &Record::declare_type_key(boost::shared_ptr<Selection> key_type) {
	ASSERT(data_->keys.size() == 0, "Declaration of TYPE key must be carried as the first.");
	data_->declare_key("TYPE", key_type, Default::obligatory(),
			"Sub-record selection.");
	return *this;
}

Record &Record::has_obligatory_type_key() {
	ASSERT( ! data_->parent_ptr_, "Record with obligatory TYPE key can't be derived.\n");
	boost::shared_ptr<Selection> sel = boost::make_shared<Selection>(type_name() + "_TYPE_selection");
	sel->add_value(0, type_name());
	sel->close();
	declare_type_key( sel );
	return *this;
}


/**********************************************************************************
 * implementation of Type::Record::RecordData
 */

Record::RecordData::RecordData(const string & type_name_in, const string & description)
:description_(description),
 type_name_(type_name_in),
 finished(false),
 closed_(false),
 derived_(false),
 auto_conversion_key_idx(-1)    // auto conversion turned off
{

}


Record::KeyIter Record::RecordData::auto_conversion_key_iter() const {
	if (auto_conversion_key_idx >= 0) return keys.begin() + auto_conversion_key_idx;
  	else return keys.end();
}



void Record::RecordData::declare_key(const string &key,
                         boost::shared_ptr<const TypeBase> type,
                         const Default &default_value, const string &description)
{
    if (finished) xprintf(PrgErr, "Declaration of key: %s in finished Record type: %s\n", key.c_str(), type_name_.c_str());

    if (key!="TYPE" && ! TypeBase::is_valid_identifier(key))
        xprintf(PrgErr, "Invalid key identifier %s in declaration of Record type: %s\n", key.c_str(), type_name_.c_str());

    KeyHash key_h = key_hash(key);
    key_to_index_const_iter it = key_to_index.find(key_h);
    if ( it == key_to_index.end() ) {
       key_to_index.insert( std::make_pair(key_h, keys.size()) );
       Key tmp_key = { (unsigned int)keys.size(), key, description, type, default_value, false};
       keys.push_back(tmp_key);
    } else {
       if (keys[it->second].derived) {
        Key tmp_key = { it->second, key, description, type, default_value, false};
        keys[ it->second ] = tmp_key;
       } else {
           xprintf(Err,"Re-declaration of the key: %s in Record type: %s\n", key.c_str(), type_name_.c_str() );
       }
    }

}



template <class KeyType>
Record &Record::declare_key(const string &key, const KeyType &type,
                        const Default &default_value, const string &description)
// this accept only lvalues - we assume that these are not local variables
{
    // ASSERT MESSAGE: The type of declared keys has to be a class derived from TypeBase.
    BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, KeyType>::value) );
    if (data_->closed_)
        xprintf(PrgErr, "Can not add key '%s' into closed record '%s'.\n", key.c_str(), type_name().c_str());

	boost::shared_ptr<const TypeBase> type_copy = boost::make_shared<KeyType>(type);
	data_->declare_key(key, type_copy, default_value, description);

    return *this;
}



template <class KeyType>
Record &Record::declare_key(const string &key, const KeyType &type,
                        const string &description)
{
    return declare_key(key,type, Default::optional(), description);
}



template <class KeyType>
AbstractRecord &AbstractRecord::declare_key(const string &key, const KeyType &type,
                        const Default &default_value, const string &description)
{
	ASSERT( false, "AbstractRecord::declare_key is not allowed!\n");
    //Record::declare_key(key, type, default_value, description);
    return *this;
}



template <class KeyType>
AbstractRecord &AbstractRecord::declare_key(const string &key, const KeyType &type,
                        const string &description)
{
	ASSERT( false, "AbstractRecord::declare_key is not allowed!\n");
    return declare_key(key,type, Default::optional(), description);
}


// explicit instantiation of template methods

#define RECORD_DECLARE_KEY(TYPE) \
template Record & Record::declare_key<TYPE>(const string &key, const TYPE &type, const Default &default_value, const string &description); \
template Record & Record::declare_key<TYPE>(const string &key, const TYPE &type, const string &description); \
template AbstractRecord &AbstractRecord::declare_key<TYPE>(const string &key, const TYPE &type, const Default &default_value, const string &description); \
template AbstractRecord &AbstractRecord::declare_key<TYPE>(const string &key, const TYPE &type, const string &description)


RECORD_DECLARE_KEY(String);
RECORD_DECLARE_KEY(Integer);
RECORD_DECLARE_KEY(Double);
RECORD_DECLARE_KEY(Bool);
RECORD_DECLARE_KEY(FileName);
RECORD_DECLARE_KEY(Selection);
RECORD_DECLARE_KEY(Array);
RECORD_DECLARE_KEY(Record);
RECORD_DECLARE_KEY(AbstractRecord);
RECORD_DECLARE_KEY(AdHocAbstractRecord);



/************************************************
 * implementation of AbstractRecord
 */

AbstractRecord::AbstractRecord()
: child_data_( boost::make_shared<ChildData>( "EmptyAbstractRecord", "" ) )
{
	attributes_ = boost::make_shared<attribute_map>();
	close();
	finish();
}



AbstractRecord::AbstractRecord(const AbstractRecord& other)
: child_data_(other.child_data_)
{
	attributes_ = other.attributes_;
}



AbstractRecord::AbstractRecord(const string & type_name_in, const string & description)
: child_data_( boost::make_shared<ChildData>( type_name_in, description ) )
{
	attributes_ = boost::make_shared<attribute_map>();
}


TypeBase::TypeHash AbstractRecord::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, "Abstract");
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, child_data_->description_);
    if (child_data_->element_input_selection != nullptr) {
    	boost::hash_combine(seed, child_data_->element_input_selection->content_hash());
    }
    //for( Record &key : child_data_->list_of_childs) {
    //    boost::hash_combine(seed, key.content_hash() );
    //}
    return seed;
}



void AbstractRecord::add_descendant(const Record &subrec)
{
    ASSERT( child_data_->closed_, "Can not add descendant to AbstractRecord that is not closed.\n");

    child_data_->selection_of_childs->add_value(child_data_->list_of_childs.size(), subrec.type_name());
    child_data_->list_of_childs.push_back(subrec);
}



AbstractRecord & AbstractRecord::allow_auto_conversion(const string &type_default) {
    if (child_data_->closed_) xprintf(PrgErr, "Can not specify default value for TYPE key as the AbstractRecord '%s' is closed.\n", type_name().c_str());
    child_data_->selection_default_=Default(type_default); // default record is closed; other constructor creates the zero item
    return *this;
}



void AbstractRecord::no_more_descendants()
{
    child_data_->selection_of_childs->close();
    // check validity of possible default value of TYPE key
    if ( !child_data_->selection_default_.is_obligatory() ) { // skip for empty records
        if ( child_data_->selection_default_.has_value_at_declaration() ) {
            try {
                child_data_->selection_of_childs->valid_default( child_data_->selection_default_.value() );
            } catch (ExcWrongDefault & e) {
                xprintf(PrgErr, "Default value '%s' for TYPE key do not match any descendant of AbstractRecord '%s'.\n", child_data_->selection_default_.value().c_str(), type_name().c_str());
            }
        }
    }
    if (! finish()) xprintf(PrgErr, "Can not finish AbstractRecord when calling no_more_descendants.\n");
}


bool AbstractRecord::valid_default(const string &str) const
{
    if ( !child_data_->selection_default_.is_obligatory() )  { // skip for empty records
        if (! child_data_->selection_of_childs->is_finished()) return false;
        if ( child_data_->selection_default_.has_value_at_declaration() ) return get_default_descendant()->valid_default(str);
    }
    THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(this->type_name()));
}


const Record  & AbstractRecord::get_descendant(const string& name) const
{
    ASSERT( child_data_->selection_of_childs->is_finished(), "Can not get descendant of unfinished AbstractType\n");
    return get_descendant( child_data_->selection_of_childs->name_to_int(name) );
}



const Record  & AbstractRecord::get_descendant(unsigned int idx) const
{
    ASSERT( child_data_->selection_of_childs->is_finished(), "Can not get descendant of unfinished AbstractType\n");
    ASSERT( idx < child_data_->list_of_childs.size() , "Size mismatch.\n");
    return child_data_->list_of_childs[idx];
}



const Record * AbstractRecord::get_default_descendant() const {
    if ( !child_data_->selection_default_.is_obligatory() )  { // skip for empty records
        if ( child_data_->selection_default_.has_value_at_declaration() ) {
            return &( get_descendant( child_data_->selection_default_.value() ) );
        }
    }
    return NULL;
}



const Selection  & AbstractRecord::get_type_selection() const
{
    return * child_data_->selection_of_childs;
}


unsigned int AbstractRecord::child_size() const {
    return child_data_->list_of_childs.size();
}


int AbstractRecord::add_child(const Record &subrec)
{
	add_descendant(subrec);
	//subrec.make_copy_keys(*this);

	return 1;
}


bool AbstractRecord::finish() {
	if (child_data_->finished_) return true;

	ASSERT(child_data_->closed_, "Finished AbstractRecord '%s' must be closed!", this->type_name().c_str());

	child_data_->finished_ = true;
	no_more_descendants();

	return (child_data_->finished_);
}


AbstractRecord &AbstractRecord::close() {
	child_data_->closed_=true;
    return *( Input::TypeRepository<AbstractRecord>::getInstance().add_type( *this ).get() );
}


AbstractRecord &AbstractRecord::set_element_input(const Selection * element_input) {
	if (element_input != NULL ) {
		child_data_->element_input_selection = element_input;
	}
	return *this;
}


bool AbstractRecord::is_finished() const {
    return child_data_->finished_;
}


bool AbstractRecord::is_closed() const {
	return child_data_->closed_;
}


string AbstractRecord::type_name() const {
    return child_data_->type_name_;
}


string AbstractRecord::full_type_name() const {
    return this->type_name();
}


Default &AbstractRecord::get_selection_default() const {
	return child_data_->selection_default_;
}



AbstractRecord::ChildDataIter AbstractRecord::begin_child_data() const {
    return child_data_->list_of_childs.begin();
}

AbstractRecord::ChildDataIter AbstractRecord::end_child_data() const {
    return child_data_->list_of_childs.end();
}


/************************************************
 * implementation of AdHocAbstractRecord
 */

AdHocAbstractRecord::AdHocAbstractRecord(const AbstractRecord &ancestor)
: AbstractRecord("Derived AdHocAbstractRecord", "This description doesn't have print out.")
{
	parent_data_ = ancestor.child_data_;
	parent_name_ = ancestor.type_name();
	tmp_ancestor_ = NULL;

	//test default descendant of ancestor
	const Record * default_desc = ancestor.get_default_descendant();
	if (default_desc) {
		allow_auto_conversion( default_desc->type_name() );
	}

	this->close();

}


AdHocAbstractRecord &AdHocAbstractRecord::add_child(const Record &subrec)
{
	AbstractRecord::add_descendant(subrec);

	return *this;
}


bool AdHocAbstractRecord::finish()
{
	if (child_data_->finished_) return true;

	if (tmp_ancestor_ != 0) {
		const_cast<AbstractRecord *>(tmp_ancestor_)->finish();

        parent_data_ = tmp_ancestor_->child_data_;
        parent_name_ = tmp_ancestor_->type_name();

		//test default descendant of ancestor
		const Record * default_desc = tmp_ancestor_->get_default_descendant();
		if (default_desc) {
			allow_auto_conversion( default_desc->type_name() );
		}

		tmp_ancestor_ = NULL;
	}

	while (unconstructed_childs_.size()) {
		const Record * rec = *(unconstructed_childs_.begin());

	    child_data_->selection_of_childs->add_value(child_data_->list_of_childs.size(), rec->type_name());
	    child_data_->list_of_childs.push_back(*rec);
	    unconstructed_childs_.pop_front();
	}

	for (AbstractRecord::ChildDataIter it = parent_data_->list_of_childs.begin(); it != parent_data_->list_of_childs.end(); ++it) {
	    child_data_->selection_of_childs->add_value(child_data_->list_of_childs.size(), (*it).type_name());
	    child_data_->list_of_childs.push_back(*it);
	}

	return AbstractRecord::finish();
}


} // closing namespace Type
} // closing namespace Input


