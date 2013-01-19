/*
 * input_type.cc
 *
 *  Created on: Mar 29, 2012
 *      Author: jb
 */


#include "input_type.hh"

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
{
    if (type_ == no_default_obligatory_type) value_="OBLIGATORY";
    if (type_ == no_default_optional_type) value_="OPTIONAL";
}


/**********************************************************************************
 * implementation of Type::Record
 */


Record::Record()
: data_( boost::make_shared<RecordData> ("EmptyRecord","") )
{
    finish();
}



Record::Record(const Record & other)
: TypeBase( other ), data_(other.data_)
{
    ASSERT( TypeBase::was_constructed(&other), "Trying to copy non-constructed Record.\n");
}



Record::Record(const string & type_name_in, const string & description)
: data_( boost::make_shared<RecordData>(type_name_in, description) )

{
    TypeBase::lazy_type_list().push_back( boost::make_shared<Record>( *this ) );
}



Record &Record::allow_auto_conversion(const string &from_key) {
    ASSERT(data_->auto_conversion_key_idx == -1, "Can not use key %s for auto conversion, the key is already set.", from_key.c_str());
    data_->auto_conversion_key_idx = 0;
    data_->auto_conversion_key=from_key;

    return *this;
}



void Record::make_derive_from(AbstractRecord &parent) const {
    if (data_->derived_) return;

    parent.finish();
    parent.add_descendant(*this);

    // copy keys form parent
    std::vector<Key>::iterator it = data_->keys.begin();
    int n_inserted = 0;
    for(KeyIter pit=parent.data_->keys.begin(); pit != parent.data_->keys.end(); ++pit) {
        Key tmp_key=*pit;    // make temporary copy of the key
        KeyHash key_h = key_hash(tmp_key.key_);

        tmp_key.derived = true;

        // we have to copy TYPE also since there should be place in storage for it
        // however we change its Default to name of actual Record
        if (tmp_key.key_=="TYPE")
            tmp_key.default_=Default( type_name() );

        // check for duplicate keys, save only the key derived by the descendant
        RecordData::key_to_index_const_iter kit = data_->key_to_index.find(key_h);
        if (kit != data_->key_to_index.end()) {
            Key *k = &(data_->keys[kit->second+n_inserted]);
            tmp_key = { tmp_key.key_index, k->key_, k->description_, k->type_, k->p_type, k->default_, false };
            k->key_ = "";
        }

        data_->key_to_index[key_h] = tmp_key.key_index;

        it = data_->keys.insert(it, tmp_key)+1;
        n_inserted++;
    }
    // delete duplicate keys and update key indices
    for (int i=0; i<data_->keys.size(); i++) {
        if (data_->keys[i].key_.compare("") == 0) {
            data_->keys.erase( data_->keys.begin()+i);
            i--;
        } else {
            data_->keys[i].key_index = i;
            data_->key_to_index[key_hash( data_->keys[i].key_)] = i;
        }
    }

    data_->derived_ = true;
}



Record &Record::derive_from(AbstractRecord &parent) {
    if (TypeBase::was_constructed(&parent)) {
        data_->parent_ptr_=boost::make_shared<AbstractRecord>(parent);
        data_->p_parent_ = NULL;
    } else { //postponed
        data_->p_parent_ = &parent;
    }

	return *this;
}



bool Record::is_finished() const {
    return data_->finished;
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
}




bool Record::finish() const
{
	if (data_->finished) return true;

	close();
    // Set correctly data_->parent_ptr
    if (data_->p_parent_ != 0 ) {
        if (TypeBase::was_constructed( data_->p_parent_))  data_->parent_ptr_=boost::make_shared<AbstractRecord>( * data_->p_parent_ );
        else return false;
        data_->p_parent_ = NULL;
    }
    // finish inheritance if parent is non-null
    if (data_->parent_ptr_) make_derive_from(* (data_->parent_ptr_));

    // Finish declare_key():
    data_->finished = true;
    for (vector<Key>::iterator it=data_->keys.begin(); it!=data_->keys.end(); it++)
    {
        // make our own copy of type object allocated at heap (could be expensive, but we don't care)
        if (it->p_type != 0) {
            if (! was_constructed(it->p_type)) return ( data_->finished=false );
            //if (! it->p_type->finish()) return false;

            if ( dynamic_cast<const AbstractRecord *>(it->p_type) != 0 ) {
                const AbstractRecord *ar = dynamic_cast<const AbstractRecord *>(it->p_type);
                it->type_ = boost::make_shared<const AbstractRecord>(*ar);
                it->p_type = 0;
            } else if (dynamic_cast<const Record *>(it->p_type) != 0) {
                const Record *r= dynamic_cast<const Record *>(it->p_type);
                it->type_ = boost::make_shared<const Record>(*r);
                it->p_type = 0;
            } else if (dynamic_cast<const Selection *>(it->p_type) != 0) {
                const Selection *s = dynamic_cast<const Selection *>(it->p_type);
                it->type_ = boost::make_shared<const Selection>(*s);
                it->p_type = 0;
            } else {
                xprintf(PrgErr,"Raw pointer to key '%s' of unknown type in record '%s'.\n", it->key_.c_str(), type_name().c_str());
            }
        }
        if (it->key_ != "TYPE") {
            data_->finished = data_->finished && it->type_->finish();

            // we check once more even keys that was already checked, otherwise we have to store
            // result of validity check in every key
            check_key_default_value(it->default_, *(it->type_), it->key_);
        }
    }

    // Check default values of autoconvertible records
    if (data_->auto_conversion_key_idx != -1 ) {
        data_->auto_conversion_key_idx=key_index(data_->auto_conversion_key);

        // check that all other keys have default values
        for(KeyIter it=data_->keys.begin(); it != data_->keys.end(); ++it) {
            if (! it->default_.has_value_at_declaration() && it->key_index != data_->auto_conversion_key_idx)
                xprintf(PrgErr, "Finishing Record auto convertible from the key '%s', but other key: '%s' has no default value.\n",
                        data_->auto_conversion_key_iter()->key_.c_str(), it->key_.c_str());
        }
    }

    return (data_->finished);
}



const Record &Record::close() const {
    data_->closed_=true;
    return *this;
}



string Record::type_name() const {
    return data_->type_name_;
}


string Record::description() const  {
    return data_->description_;
}



bool Record::valid_default(const string &str) const
{
    if (data_->auto_conversion_key_idx >=0) {
        unsigned int idx=key_index(data_->auto_conversion_key);
        if ( data_->keys[idx].type_ ) return data_->keys[idx].type_->valid_default(str);
        else return false;
    } else {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(this->type_name()));
    }
}

/*
void  Record::reset_doc_flags() const {
    data_->reset_doc_flags();
}
*/


bool Record::operator==(const TypeBase &other) const
{ return  typeid(*this) == typeid(other) &&
                 (type_name() == static_cast<const Record *>(&other)->type_name() );
}


Record::KeyIter Record::auto_conversion_key_iter() const {
    finished_check();
    return data_->auto_conversion_key_iter();
}


/**********************************************************************************
 * implementation of Type::Record::RecordData
 */

Record::RecordData::RecordData(const string & type_name_in, const string & description)
:description_(description),
 type_name_(type_name_in),
 p_parent_(0),
 made_extensive_doc(false),
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


/*
void  Record::RecordData::reset_doc_flags() const {
    made_extensive_doc=false;
    for(KeyIter it = keys.begin(); it!=keys.end(); ++it) {
        it->type_->reset_doc_flags();
    }
}
*/



void Record::RecordData::declare_key(const string &key,
                         boost::shared_ptr<const TypeBase> type,
                         const TypeBase *type_temporary,
                         const Default &default_value, const string &description)
{
    if (finished) xprintf(PrgErr, "Declaration of key: %s in finished Record type: %s\n", key.c_str(), type_name_.c_str());

    if (key!="TYPE" && ! TypeBase::is_valid_identifier(key))
        xprintf(PrgErr, "Invalid key identifier %s in declaration of Record type: %s\n", key.c_str(), type_name_.c_str());

    KeyHash key_h = key_hash(key);
    key_to_index_const_iter it = key_to_index.find(key_h);
    if ( it == key_to_index.end() ) {
       key_to_index.insert( std::make_pair(key_h, keys.size()) );
       Key tmp_key = { (unsigned int)keys.size(), key, description, type, type_temporary, default_value, false};
       keys.push_back(tmp_key);
    } else {
       if (keys[it->second].derived) {
    	   Key tmp_key = { it->second, key, description, type, type_temporary, default_value, false};
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
    if ( (boost::is_base_of<Record, KeyType>::value ||
          boost::is_base_of<Selection, KeyType>::value)
         && ! TypeBase::was_constructed(&type) ) {

        data_->declare_key(key, boost::shared_ptr<const TypeBase>(), &type, default_value, description);
    } else {
        // for Array, Double, Integer, we assume no static variables
        if (type.finish()) check_key_default_value(default_value, type, key);

        boost::shared_ptr<const TypeBase> type_copy = boost::make_shared<KeyType>(type);
        data_->declare_key(key, type_copy, NULL, default_value, description);
    }

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
    Record::declare_key(key, type, default_value, description);
    return *this;
}



template <class KeyType>
AbstractRecord &AbstractRecord::declare_key(const string &key, const KeyType &type,
                        const string &description)
{
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



/************************************************
 * implementation of AbstractRecord
 */

AbstractRecord::AbstractRecord()
: Record(), child_data_( boost::make_shared<ChildData>( "EmptyAbstractRecord_TYPE_selection" ) )
{
}



AbstractRecord::AbstractRecord(const AbstractRecord& other)
: Record(other), child_data_(other.child_data_)
{
    ASSERT( TypeBase::was_constructed(&other), "Trying to copy non-constructed Record.\n");
}



AbstractRecord::AbstractRecord(const string & type_name_in, const string & description)
: Record(type_name_in, description),
  child_data_( boost::make_shared<ChildData>( type_name_in + "_TYPE_selection" ) )
{
    // declare very first item of any descendant
    data_->declare_key("TYPE", boost::shared_ptr<Selection>(), child_data_->selection_of_childs.get(), Default::obligatory(),
                 "Sub-record selection.");

    TypeBase::lazy_type_list().push_back( boost::make_shared<AbstractRecord>( *this ) );
}



void AbstractRecord::add_descendant(const Record &subrec)
{
    ASSERT( data_->closed_, "Can not add descendant to AbstractRecord that is not closed.\n");

    child_data_->selection_of_childs->add_value(child_data_->list_of_childs.size(), subrec.type_name());
    child_data_->list_of_childs.push_back(subrec);
}



AbstractRecord & AbstractRecord::allow_auto_conversion(const string &type_default) {
    ASSERT( ! data_->closed_, "Can not specify default value for TYPE key as the AbstractRecord '%s' is closed.\n", type_name().c_str());
    data_->keys[0].default_=Default(type_default); // default record is closed; other constructor creates the zero item
    return *this;
}



void AbstractRecord::no_more_descendants()
{
    child_data_->selection_of_childs->close();
    // check validity of possible default value of TYPE key
    if (data_->keys.size() > 0 ) { // skip for empty records
        Default &dflt = data_->keys[0].default_;
        if ( dflt.has_value_at_declaration() ) {
            try {
                child_data_->selection_of_childs->valid_default( dflt.value() );
            } catch (ExcWrongDefault & e) {
                xprintf(PrgErr, "Default value '%s' for TYPE key do not match any descendant of AbstractRecord '%s'.\n", data_->keys[0].default_.value().c_str(), type_name().c_str());
            }
        }
    }
    if (! finish()) xprintf(PrgErr, "Can not finish AbstractRecord when calling no_more_descendants.\n");
}


/*
void  AbstractRecord::reset_doc_flags() const {
        data_->reset_doc_flags();
        for(vector< Record >::const_iterator it=child_data_->list_of_childs.begin();
                    it!= child_data_->list_of_childs.end(); ++it)
            it->reset_doc_flags();
}
*/


const Record  & AbstractRecord::get_descendant(const string& name) const
{
    ASSERT( is_finished(), "Can not get descendant of unfinished AbstractType\n");
    return get_descendant( child_data_->selection_of_childs->name_to_int(name) );
}



const Record  & AbstractRecord::get_descendant(unsigned int idx) const
{

    ASSERT( idx < child_data_->list_of_childs.size() , "Size mismatch.\n");
    return child_data_->list_of_childs[idx];
}



const Record * AbstractRecord::get_default_descendant() const {
    if (data_->keys.size() != 0 )  { // skip for empty records
        Default &dflt = data_->keys[0].default_;
        if ( dflt.has_value_at_declaration() ) {
            return &( get_descendant( child_data_->selection_of_childs->name_to_int( dflt.value() ) ) );
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


AbstractRecord::ChildDataIter AbstractRecord::begin_child_data() const {
    child_data_->list_of_childs.begin();
}

AbstractRecord::ChildDataIter AbstractRecord::end_child_data() const {
    child_data_->list_of_childs.end();
}




} // closing namespace Type
} // closing namespace Input


