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

#include "system.hh"

#include <boost/type_traits.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>

namespace Input {
namespace Type {

using namespace std;


/**********************************************************************************
 * implementation of Type::Record
 */

Record::Record(const string & type_name_in, const string & description)
: data_( boost::make_shared<RecordData>(type_name_in, description) )

{
    finished=false;
}


void Record::derive_from(AbstractRecord parent) {

    if (data_.use_count() == 0)
            xprintf(PrgErr, "Can not inherit to empty Record handle.\n");
    if (data_->keys.size() != 0)
            xprintf(PrgErr, "Can not inherit into Record `%s`, it already has some keys declared\n.", type_name().c_str());

    parent.add_descendant(*this);
    // semi-deep copy
    for(Record::KeyIter it=parent.begin(); it != parent.end(); ++it) {
        Key tmp_key=*it;
        KeyHash key_h = key_hash(tmp_key.key_);

        data_->key_to_index.insert( std::make_pair(key_h, data_->keys.size()) );
        tmp_key.key_index=data_->keys.size();
        tmp_key.derived = true;
        data_->keys.push_back(tmp_key);
    }
}

string Record::type_name() const
{
    if (data_.use_count() == 0) return "empty_handle";
    else return data_->type_name_;
}

void  Record::reset_doc_flags() const {
    if (data_.use_count() != 0) data_->reset_doc_flags();
}

std::ostream& Record::documentation(std::ostream& stream, bool extensive, unsigned int pad) const
{
    if (! finished) xprintf(PrgErr, "Can not provide documentation of unfinished Record type: %s\n", type_name().c_str());

    return data_->documentation(stream, extensive, pad);
}

/**********************************************************************************
 * implementation of Type::Record::RecordData
 */

Record::RecordData::RecordData(const string & type_name_in, const string & description)
:description_(description),
 type_name_(type_name_in),
 made_extensive_doc(false)
{}

std::ostream& Record::RecordData::documentation(std::ostream& stream, bool extensive, unsigned int pad) const
{

    if (! extensive) {

        // Short description
        stream << "Record '" << type_name_ << "' with "<< keys.size() << " keys";
    } else if ( ! made_extensive_doc) {

        // Extensive description
        made_extensive_doc=true;

        // header
        stream << endl;
        stream << ""
               << "Record '" << type_name_ << "' with "<< keys.size() << " keys.";
        write_description(stream, description_, pad);
        stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
        // keys
        for(KeyIter it = keys.begin(); it!=keys.end(); ++it) {
            stream << setw(pad + 4) << ""
                   << it->key_ << " = <" << it->default_.value() << "> is ";
            it->type_->documentation( stream , false, pad +4 ); // short description of the type of the key
            write_description(stream, it->description_, pad + 4); // description of the key on further lines

        }
        stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ')
        << " " << type_name_ << endl;

        // Full documentation of embedded record types.
        for(KeyIter it = keys.begin(); it!=keys.end(); ++it) {
            it->type_->documentation(stream, true, 0);
        }

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
                         const DefaultValue &default_value, const string &description)
{
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

/************************************************
 * implementation of AbstractRecord
 */

AbstractRecord::AbstractRecord(const string & type_name_in, const string & description)
: Record(type_name_in, description),
  child_data_( boost::make_shared<ChildData>( type_name_in + "_selection" ) )
{
    // make our own copy of type object allocated at heap (could be expensive, but we don't care)
    boost::shared_ptr< Selection<unsigned int> > type_copy = boost::make_shared< Selection<unsigned int> >(child_data_->selection_of_childs);

//    data_->declare_key("TYPE", type_copy, DefaultValue(DefaultValue::obligatory),
//                 "Sub-record selection.");

    finished=false;
}



void AbstractRecord::add_descendant(const Record &subrec)
{
    if (!finished)
        xprintf(PrgErr, "Can not add descendant to unfinished AbstractType.\n");

    child_data_->selection_of_childs.add_value(child_data_->list_of_childs.size(), subrec.type_name());
    child_data_->list_of_childs.push_back(subrec);
}



void AbstractRecord::no_more_descendants()
{
    if (child_data_.use_count() == 0)
            xprintf(PrgErr, "Can not close empty AbstractRecord handle.\n");
    child_data_->selection_of_childs.finish();

}



void  AbstractRecord::reset_doc_flags() const {
    if (data_.use_count() != 0) {
        data_->reset_doc_flags();
        for(vector< Record >::const_iterator it=child_data_->list_of_childs.begin();
                    it!= child_data_->list_of_childs.end(); ++it)
            it->reset_doc_flags();
    }
}



std::ostream& AbstractRecord::documentation(std::ostream& stream, bool extensive, unsigned int pad) const
{
    if (! finished) xprintf(PrgErr, "Can not provide documentation of unfinished Record type: %s\n", type_name().c_str());

    if (! extensive) {

        // Short description
        stream << "AbstractRecord '" << type_name() << "' with "<< child_data_->list_of_childs.size() << " descendants.";
    } else if ( ! data_->made_extensive_doc) {

        // Extensive description
        data_->made_extensive_doc=true;

        // header
        stream << endl;
        stream << ""
               << "AbstractRecord '" << type_name() << "' with "<< child_data_->list_of_childs.size() << " descendants.";
        write_description(stream, data_->description_, pad);
        stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
        // descendants
        for(vector< Record >::const_iterator it=child_data_->list_of_childs.begin();
                    it!= child_data_->list_of_childs.end(); ++it) {
            it->documentation( stream , false, pad +4 ); // short description of the type of the key
        }
        stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ')
        << " " << type_name() << endl;

        // Full documentation of embedded record types.
        for(vector< Record >::const_iterator it=child_data_->list_of_childs.begin();
                    it!= child_data_->list_of_childs.end(); ++it) {
            it->documentation(stream, true, 0);
        }

    }

    return stream;
}




const Record  & AbstractRecord::get_descendant(const string& name) const
{
    unsigned int idx;

    ASSERT( finished, "Can not get descendant of unfinished AbstractType\n");
    idx = child_data_->selection_of_childs.name_to_int(name);
    ASSERT( idx < child_data_->list_of_childs.size() , "Size mismatch.\n");
    return child_data_->list_of_childs[idx];
}



} // closing namespace Type
} // closing namespace Input


