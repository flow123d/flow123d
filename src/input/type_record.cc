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


Record::Record() {
    TypeBase::insert_lazy_object(this);
}



Record::Record(const Record & other)
: TypeBase( other )
{
    ASSERT( TypeBase::was_constructed(&other), "Trying to copy non-constructed Record.\n");
    TypeBase::insert_lazy_object(this);
    data_ = other.data_;
}



Record::Record(const string & type_name_in, const string & description)
: data_( boost::make_shared<RecordData>(type_name_in, description) )

{
    TypeBase::insert_lazy_object(this);
    TypeBase::lazy_type_list().push_back( boost::make_shared<Record>( *this ) );
}



Record &Record::allow_auto_conversion(const string &from_key) {
    empty_check();
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
        // however we change its Default to optional()
        if (tmp_key.key_=="TYPE")
            tmp_key.default_=Default::optional();

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

    empty_check();
    if (TypeBase::was_constructed(&parent)) {
        data_->parent_ptr_=boost::make_shared<AbstractRecord>(parent);
        data_->p_parent_ = NULL;
    } else { //postponed
        data_->p_parent_ = &parent;
    }

	return *this;
}



bool Record::is_finished() const {
    empty_check();
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
	empty_check();
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
        }
        // we check once more even keys that was already checked, otherwise we have to store
        // result of validity check in every key
        check_key_default_value(it->default_, *(it->type_), it->key_);
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
    empty_check();
    data_->closed_=true;
    return *this;
}



std::ostream& Record::documentation(std::ostream& stream,DocType extensive, unsigned int pad) const
{
    if (! is_finished()) xprintf(Warn, "Printing documentation of unfinished Input::Type::Record!\n");
    return data_->documentation(stream, extensive, pad);
}



string Record::type_name() const
{
    if (data_.use_count() == 0) return "empty_handle";
    else return data_->type_name_;
}



string Record::description() const
{
    if (data_.use_count() == 0) return "empty_handle";
    else return data_->description_;
}



bool Record::valid_default(const string &str) const
{
    empty_check();
    if (data_->auto_conversion_key_idx >=0) {
        unsigned int idx=key_index(data_->auto_conversion_key);
        if ( data_->keys[idx].type_ ) return data_->keys[idx].type_->valid_default(str);
        else return false;
    } else {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(this->type_name()));
    }
}


void  Record::reset_doc_flags() const {
    if (data_.use_count() != 0) data_->reset_doc_flags();
}



bool Record::operator==(const TypeBase &other) const
{ return  typeid(*this) == typeid(other) &&
                 (type_name() == static_cast<const Record *>(&other)->type_name() );
}


Record::KeyIter Record::auto_conversion_key_iter() const {
    empty_check();
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







std::ostream& Record::RecordData::documentation(std::ostream& stream,DocType extensive, unsigned int pad) const
{

    switch (extensive) {
    case record_key:
        // Short description
        stream << "Record '" << type_name_ << "' with "<< keys.size() << " keys";
        break;
    case full_after_record:
        // Detailed with recursion
        if (!made_extensive_doc) {

            // Extensive description
            made_extensive_doc = true;
            pad=0;

            // header
            stream << endl;
            stream << "" << "Record '" << type_name_ << "' with " << keys.size() << " keys.";
            write_description(stream, description_, pad);
            stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
            // keys
            for (KeyIter it = keys.begin(); it != keys.end(); ++it) {
                stream << setw(pad + 4) << "" << it->key_ << " = <" << it->default_.value() << "> is ";
                it->type_->documentation(stream, record_key, pad + 4); // short description of the type of the key
                write_description(stream, it->description_, pad + 4); // description of the key on further lines

            }
            stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type_name_ << endl;

            // Full documentation of embedded record types.
            for (KeyIter it = keys.begin(); it != keys.end(); ++it) {
                it->type_->documentation(stream, full_after_record, 0);
            }

        }
        break;
    case full_along:
        pad=0;

        // header
        stream << endl;
        stream << "" << "Record '" << type_name_ << "' with " << keys.size() << " keys.";
        write_description(stream, description_, pad);
        stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
        // keys
        for (KeyIter it = keys.begin(); it != keys.end(); ++it) {
            stream << setw(pad + 4) << "" << it->key_ << " = <" << it->default_.value() << "> is ";
            it->type_->documentation(stream, record_key, pad + 4); // short description of the type of the key
            write_description(stream, it->description_, pad + 4); // description of the key on further lines

        }
        stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type_name_ << endl;
        break;
    }

    return stream;
}



void  Record::RecordData::reset_doc_flags() const {
    made_extensive_doc=false;
    for(KeyIter it = keys.begin(); it!=keys.end(); ++it) {
        it->type_->reset_doc_flags();
    }
}




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



/************************************************
 * implementation of AbstractRecord
 */

AbstractRecord::AbstractRecord()
: Record(), child_data_()
{
    TypeBase::insert_lazy_object(this);
}



AbstractRecord::AbstractRecord(const AbstractRecord& other)
: Record(other)
{
    ASSERT( TypeBase::was_constructed(&other), "Trying to copy non-constructed Record.\n");

    TypeBase::insert_lazy_object(this);
    child_data_ = other.child_data_;
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
    empty_check();
    ASSERT( data_->closed_, "Can not add descendant to AbstractRecord that is not closed.\n");

    child_data_->selection_of_childs->add_value(child_data_->list_of_childs.size(), subrec.type_name());
    child_data_->list_of_childs.push_back(subrec);
}



void AbstractRecord::no_more_descendants()
{
    if (child_data_.use_count() == 0)
            xprintf(PrgErr, "Can not close empty AbstractRecord handle.\n");
    child_data_->selection_of_childs->close();
    empty_check();
    if (! finish()) xprintf(PrgErr, "Can not finish AbstractRecord when calling no_more_descendants.\n");
}



void  AbstractRecord::reset_doc_flags() const {
    if (data_.use_count() != 0) {
        data_->reset_doc_flags();
        for(vector< Record >::const_iterator it=child_data_->list_of_childs.begin();
                    it!= child_data_->list_of_childs.end(); ++it)
            it->reset_doc_flags();
    }
}



std::ostream& AbstractRecord::documentation(std::ostream& stream,DocType extensive, unsigned int pad) const
{
    if (! is_finished()) xprintf(Warn, "Printing documentation of unfinished Input::Type:AbstractRecord!\n");

    switch (extensive) {
    case record_key:
        // Short description
        stream << "AbstractRecord '" << type_name() << "' with "<< child_data_->list_of_childs.size() << " descendants.";
        break;
    case full_after_record:
        if (!data_->made_extensive_doc) {

            // Extensive description
            data_->made_extensive_doc = true;
            pad=0;

            // header
            stream << endl;
            stream << "" << "AbstractRecord '" << type_name() << "' with " << child_data_->list_of_childs.size() << " descendants.";
            write_description(stream, data_->description_, pad);
            stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
            // descendants
            for (vector<Record>::const_iterator it = child_data_->list_of_childs.begin(); it != child_data_->list_of_childs.end();
                    ++it) {
                stream << setw(pad + 4) << "";
                it->documentation(stream, record_key, 0); // short description of the type of the key
                stream << endl;
            }
            stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type_name() << endl;

            // Full documentation of embedded record types.
            for (vector<Record>::const_iterator it = child_data_->list_of_childs.begin(); it != child_data_->list_of_childs.end();
                    ++it) {
                it->documentation(stream, full_after_record, 0);
            }

        }
        break;
    case full_along:
        // Extensive description
        data_->made_extensive_doc = true;
        pad=0;

        // header
        stream << endl;
        stream << "" << "AbstractRecord '" << type_name() << "' with " << child_data_->list_of_childs.size() << " descendants.";
        write_description(stream, data_->description_, pad);
        stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
        // descendants
        for (vector<Record>::const_iterator it = child_data_->list_of_childs.begin(); it != child_data_->list_of_childs.end();
                ++it) {
            stream << setw(pad + 4) << "";
            it->documentation(stream, record_key, 0); // short description of the type of the key
            stream << endl;
        }
        stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type_name() << endl;
        break;
    }
    return stream;
}




const Record  & AbstractRecord::get_descendant(const string& name) const
{
    ASSERT(child_data_.use_count() != 0, "Wrong use of an empty AbstractRecord.");
    ASSERT( is_finished(), "Can not get descendant of unfinished AbstractType\n");
    return get_descendant( child_data_->selection_of_childs->name_to_int(name) );
}



const Record  & AbstractRecord::get_descendant(unsigned int idx) const
{
    ASSERT(child_data_.use_count() != 0, "Wrong use of an empty AbstractRecord.");

    ASSERT( idx < child_data_->list_of_childs.size() , "Size mismatch.\n");
    return child_data_->list_of_childs[idx];
}


const Selection  & AbstractRecord::get_type_selection() const
{
    ASSERT(child_data_.use_count() != 0, "Wrong use of an empty AbstractRecord.");
    return * child_data_->selection_of_childs;
}



} // closing namespace Type
} // closing namespace Input


