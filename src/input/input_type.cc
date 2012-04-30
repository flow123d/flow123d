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

std::ostream& Record::documentation(std::ostream& stream, bool extensive, unsigned int pad) const
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

void  Record::reset_doc_flags() const {
    made_extensive_doc=false;
    for(KeyIter it = keys.begin(); it!=keys.end(); ++it) {
        it->type_->reset_doc_flags();
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

std::ostream& Integer::documentation(std::ostream& stream, bool extensive, unsigned int pad)  const {
    if (extensive) return stream;
    stream << "Integer in [" << lower_bound_ << ", " << upper_bound_ << "]";
    return stream;
}

std::ostream& Double::documentation(std::ostream& stream, bool extensive, unsigned int pad)  const {
    if (extensive) return stream;
    stream << "Double in [" << lower_bound_ << ", " << upper_bound_ << "]";
    return stream;
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





std::ostream& String::documentation(std::ostream& stream, bool extensive, unsigned int pad) const {

    if (extensive) return stream;
    stream << "String (generic)";
    return stream;
}






} // closing namespace Type
} // closing namespace Input


