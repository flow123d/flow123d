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


/*******************************************************************
 * implementation of DefaultValue
 */

DefaultValue::DefaultValue()
: value_(), type_(optional)
{}

DefaultValue::DefaultValue(const std::string & value)
: value_(value), type_(declaration)
{}

DefaultValue::DefaultValue(enum DefaultType type)
: value_(), type_(type)
{
    if (type == declaration) {
        xprintf(Err, "Can not construct DefaultValue with type 'declaration' without providing the default value.\n");
    }
}


/*******************************************************************
 * implementation of TypeBase
 */

std::ostream& TypeBase::write_description(std::ostream& stream, const string& str, unsigned int pad) {
    boost::tokenizer<boost::char_separator<char> > line_tokenizer(str, boost::char_separator<char>("\n"));
    boost::tokenizer<boost::char_separator<char> >::iterator tok;

    // Up to first \n without padding.
    stream << endl;

        // For every \n add padding at beginning of the nex line.
        for(tok = line_tokenizer.begin(); tok != line_tokenizer.end(); ++tok) {
            stream << setw(pad) << "" << "# "
                    << *tok << endl;
        }
    return stream;
}

bool TypeBase::is_valid_identifier(const string& key) {
  namespace ba = boost::algorithm;
  return ba::all( key, ba::is_lower() || ba::is_digit() || ba::is_any_of("_") );
}


std::ostream& operator<<(std::ostream& stream, const TypeBase& type) {
    return type.documentation(stream);
}


/**********************************************************************************
 * implementation of Type::Record
 */

Record::Record(const string & type_name, const string & description)
: data_( boost::make_shared<RecordData>(description, type_name) )

{
    finished=false;
}

Record::Record( AbstractRecord parent, const string & type_name, const string & description)
: data_( boost::make_shared<RecordData>(description, type_name) )
{
    finished=false;
    parent.add_descendant(*this);
    // semi-deep copy
    for(Record::KeyIter it=parent.begin(); it != parent.end(); ++it) {
        //...
    }

}

const string &Record::type_name() const
{ return data_->type_name_;}

void  Record::reset_doc_flags() const {
    if (data_.use_count() == 0) xprintf(PrgErr, "Can not reset flags of an empty Record proxy.\n");
    data_->reset_doc_flags();
}

std::ostream& Record::documentation(std::ostream& stream, bool extensive, unsigned int pad) const
{
    if (! finished) xprintf(PrgErr, "Can not provide documentation of unfinished Record type: %s\n", type_name().c_str());
    if (data_.use_count() == 0) xprintf(PrgErr, "Can not document an empty Record proxy.\n");

    return data_->documentation(stream, extensive, pad);
}

/**********************************************************************************
 * implementation of Type::Record::RecordData
 */

Record::RecordData::RecordData(const string & type_name, const string & description)
:description_(description),
 type_name_(type_name),
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
       Key tmp_key = { it->second, key, description, type, default_value, false};
       keys[ it->second ] = tmp_key;

       xprintf(Err,"Re-declaration of the key: %s in Record type: %s\n", key.c_str(), type_name_.c_str() );
    }
}


/**********************************************************************************
 * implementation of Type::Array
 */


std::ostream& Array::documentation(std::ostream& stream, bool extensive, unsigned int pad) const {
    if (extensive) {
        type_of_values_->documentation(stream, true, pad);
    } else {
        stream << "Array, size limits: [" << lower_bound_ << ", " << upper_bound_ << "] of type: " << endl;
        stream << setw(pad+4) << "";
        type_of_values_->documentation(stream, false, pad+4);
    }
    return stream;
}

void  Array::reset_doc_flags() const {
    type_of_values_->reset_doc_flags();
}

const string &Array::type_name() const
{ return "array_of_" + type_of_values_->type_name();}

/**********************************************************************************
 * implementation of Type::Scalar ... and descendants.
 */
void  Scalar::reset_doc_flags() const
{}

std::ostream& Bool::documentation(std::ostream& stream, bool extensive, unsigned int pad)  const {
    if (extensive) return stream;
    stream << "Bool";
    return stream;
}

const string &Bool::type_name() const {
    return "Bool";
}

std::ostream& Integer::documentation(std::ostream& stream, bool extensive, unsigned int pad)  const {
    if (extensive) return stream;
    stream << "Integer in [" << lower_bound_ << ", " << upper_bound_ << "]";
    return stream;
}

const  string &Integer::type_name() const {
    return "Integer";
}


std::ostream& Double::documentation(std::ostream& stream, bool extensive, unsigned int pad)  const {
    if (extensive) return stream;
    stream << "Double in [" << lower_bound_ << ", " << upper_bound_ << "]";
    return stream;
}

const string &Double::type_name() const {
    return "Double";
}


std::ostream& FileName::documentation(std::ostream& stream, bool extensive, unsigned int pad)  const {
    if (extensive) return stream;

    stream << "FileName of ";
    switch (type_) {
    case input_file:
        stream << "input file";
        break;
    case output_file:
        stream << "output file";
        break;
    default:
        stream << "file with unknown type";
        break;
    }
    return stream;
}

const string &FileName::type_name() const {
    switch (type_) {
    case input_file:
        return "FileName_input";
    case output_file:
        return "FileName_output";
    default:
        return "FileName";
    }
}




std::ostream& String::documentation(std::ostream& stream, bool extensive, unsigned int pad) const {

    if (extensive) return stream;
    stream << "String (generic)";
    return stream;
}

const string &String::type_name() const {
    return "String";
}






} // closing namespace Type
} // closing namespace Input


