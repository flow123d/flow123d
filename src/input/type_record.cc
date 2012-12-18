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

Record::Record(const string & type_name_in, const string & description)
: data_( boost::make_shared<RecordData>(type_name_in, description) )

{}


Record::Record(boost::shared_ptr<RecordData> data_ptr)
: data_(data_ptr)
{}



Record &Record::allow_auto_conversion(const string &from_key) {
    empty_check();
    ASSERT(data_->auto_conversion_key == -1, "Can not use key %s for auto conversion, the key is already set.", from_key.c_str());
    data_->auto_conversion_key = key_index(from_key);

    return *this;
}


Record &Record::derive_from(AbstractRecord &parent) {

    empty_check();
    if (data_->keys.size() != 0)
            THROW( ExcDeriveNonEmpty() << EI_RecordName(data_->parent_->type_name()) << EI_Record(*this) );

    // Reserve the first key for the TYPE selection.
    // The reference to the child list will be updated in finish().
    data_->declare_key_reference("TYPE", 0, Default::obligatory(),
                 "Sub-record selection.");

	data_->parent_ = &parent;
	data_->descendant_data_ = data_;

	return *this;
}





bool Record::is_finished() const {
    empty_check();
    return data_->finished;
}


void Record::finish()
{
	empty_check();
	data_->finish();
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



void Record::valid_default(const string &str) const
{
    empty_check();
    if (data_->auto_conversion_key >=0) {
        this->auto_conversion_key_iter()->type_->valid_default(str);
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
    return data_->auto_conversion_key_iter();
}


/**********************************************************************************
 * implementation of Type::Record::RecordData
 */

Record::RecordData::RecordData(const string & type_name_in, const string & description)
:description_(description),
 type_name_(type_name_in),
 parent_(0),
 made_extensive_doc(false),
 finished(false),
 auto_conversion_key(-1)    // auto conversion turned off
{
	LazyTypes::instance().addType(this);
}


Record::KeyIter Record::RecordData::auto_conversion_key_iter() const {
	if (auto_conversion_key >= 0) return keys.begin() + auto_conversion_key;
  	else return keys.end();
}



void Record::RecordData::finish() {

    if (finished) return;

    // Finish derive_from():
    if (parent_ != 0)
    {
    	parent_->finish();

    	Record tmp_rec(descendant_data_);
    	parent_->add_descendant(tmp_rec);
		// copy keys form parent, we have to copy TYPE also since there should be place in storage for it
		// however we change its Default to optional()
		for(KeyIter it=parent_->begin(); it != parent_->end(); ++it) {
			Key tmp_key=*it;    // make temporary copy of the key
			KeyHash key_h = key_hash(tmp_key.key_);

			tmp_key.derived = true;

			// the TYPE key is placed to the beginning of the vector keys,
			// while the other keys are appended to the end.
			if (tmp_key.key_=="TYPE") {
				tmp_key.default_=Default::optional();
				key_to_index.insert( std::make_pair(key_h, 0) );
				tmp_key.key_index=0;
				keys[0] = tmp_key;
			} else {
				key_to_index.insert( std::make_pair(key_h, keys.size()) );
				tmp_key.key_index=keys.size();
				keys.push_back(tmp_key);
			}
		}

		parent_ = 0;
    }



    // Finish declare_key():
    for (vector<Key>::iterator it=keys.begin(); it!=keys.end(); it++)
    {
		// make our own copy of type object allocated at heap (could be expensive, but we don't care)
		if (it->p_type != 0) {
			if (dynamic_cast<const AbstractRecord *>(it->p_type) != 0) {
				it->type_ = boost::make_shared<const AbstractRecord>(*dynamic_cast<const AbstractRecord *>(it->p_type));
				it->p_type = 0;
			} else if (dynamic_cast<const Record *>(it->p_type) != 0) {
				it->type_ = boost::make_shared<const Record>(*dynamic_cast<const Record *>(it->p_type));
				it->p_type = 0;
			} else if (dynamic_cast<const Selection *>(it->p_type) != 0) {
				it->type_ = boost::make_shared<const Selection>(*dynamic_cast<const Selection *>(it->p_type));
				it->p_type = 0;
			}
		} else if (dynamic_cast<const Array *>(it->type_.get()) != 0) {
			// Arrays may be of type Record, AbstractRecord or Selection,
			// in which case they must be finished.
			((Array *)dynamic_cast<const Array *>(it->type_.get()))->finish();
		}

		// check validity of possibly given default value
		if ( it->default_.has_value_at_declaration() ) {

			try {
				it->type_->valid_default( it->default_.value() );
			} catch (ExcWrongDefault & e) {
				e << EI_KeyName(it->key_);
				throw;
			}
		}
    }

    // Check default values
    if (auto_conversion_key != -1 ) {
        // check that all other keys have default values
        for(KeyIter it=keys.begin(); it != keys.end(); ++it) {
            if (! it->default_.has_value_at_declaration() && it->key_index != auto_conversion_key)
                xprintf(PrgErr, "Finishing Record auto convertible from the key '%s', but other key: '%s' has no default value.\n",
                        auto_conversion_key_iter()->key_.c_str(), it->key_.c_str());
        }
    }

    finished = true;
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
                         const Default &default_value, const string &description)
{
    KeyHash key_h = key_hash(key);
    key_to_index_const_iter it = key_to_index.find(key_h);
    if ( it == key_to_index.end() ) {
       key_to_index.insert( std::make_pair(key_h, keys.size()) );
       Key tmp_key = { (unsigned int)keys.size(), key, description, type, 0, default_value, false};
       keys.push_back(tmp_key);
    } else {
       if (keys[it->second].derived) {
    	   Key tmp_key = { it->second, key, description, type, 0, default_value, false};
           keys[ it->second ] = tmp_key;
       } else {
           xprintf(Err,"Re-declaration of the key: %s in Record type: %s\n", key.c_str(), type_name_.c_str() );
       }
    }
}


void Record::RecordData::declare_key_reference(const string &key,
                         const TypeBase *type,
                         const Default &default_value, const string &description)
{
    KeyHash key_h = key_hash(key);
    key_to_index_const_iter it = key_to_index.find(key_h);
    if ( it == key_to_index.end() ) {
       key_to_index.insert( std::make_pair(key_h, keys.size()) );
       Key tmp_key = { (unsigned int)keys.size(), key, description, boost::shared_ptr<const TypeBase>(), type, default_value, false };
       keys.push_back(tmp_key);
    } else {
       if (keys[it->second].derived) {
           Key tmp_key = { it->second, key, description, boost::shared_ptr<const TypeBase>(), type, default_value, false };
           keys[ it->second ] = tmp_key;
       } else {
           xprintf(Err,"Re-declaration of the key: %s in Record type: %s\n", key.c_str(), type_name_.c_str() );
       }
    }
}




/************************************************
 * implementation of AbstractRecord
 */

AbstractRecord::AbstractRecord(const string & type_name_in, const string & description)
: Record(type_name_in, description),
  child_data_( boost::make_shared<ChildData>( type_name_in + "_TYPE_selection" ) )
{

    // declare very first item of any descendant
    data_->declare_key_reference("TYPE", child_data_->selection_of_childs.get(), Default::obligatory(),
                 "Sub-record selection.");

}



void AbstractRecord::add_descendant(const Record &subrec)
{
    ASSERT( is_finished(), "Can not add descendant to unfinished AbstractType.\n");

    child_data_->selection_of_childs->add_value(child_data_->list_of_childs.size(), subrec.type_name());
    child_data_->list_of_childs.push_back(subrec);
}



void AbstractRecord::no_more_descendants()
{
    if (child_data_.use_count() == 0)
            xprintf(PrgErr, "Can not close empty AbstractRecord handle.\n");
    child_data_->selection_of_childs->finish();

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


