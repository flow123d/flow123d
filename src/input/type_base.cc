/*
 * input_type.cc
 *
 *  Created on: Mar 29, 2012
 *      Author: jb
 */



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

#include "type_base.hh"


namespace Input {
namespace Type {

using namespace std;


/*******************************************************************
 * implementation of Default
 */

Default::Default()
: value_(), type_(optional_type)
{}

Default::Default(const std::string & value)
: value_(value), type_(declaration)
{}

Default::Default(enum DefaultType type)
: value_(), type_(type)
{}


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
    type.reset_doc_flags();
    return type.documentation(stream, true,0);
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

string Array::type_name() const
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

string Bool::type_name() const {
    return "Bool";
}

std::ostream& Integer::documentation(std::ostream& stream, bool extensive, unsigned int pad)  const {
    if (extensive) return stream;
    stream << "Integer in [" << lower_bound_ << ", " << upper_bound_ << "]";
    return stream;
}

string Integer::type_name() const {
    return "Integer";
}


std::ostream& Double::documentation(std::ostream& stream, bool extensive, unsigned int pad)  const {
    if (extensive) return stream;
    stream << "Double in [" << lower_bound_ << ", " << upper_bound_ << "]";
    return stream;
}

string Double::type_name() const {
    return "Double";
}


std::ostream& FileName::documentation(std::ostream& stream, bool extensive, unsigned int pad)  const {
    if (extensive) return stream;

    stream << "FileName of ";
    switch (type_) {
    case ::FilePath::input_file:
        stream << "input file";
        break;
    case ::FilePath::output_file:
        stream << "output file";
        break;
    default:
        stream << "file with unknown type";
        break;
    }
    return stream;
}

string FileName::type_name() const {
    switch (type_) {
    case ::FilePath::input_file:
        return "FileName_input";
    case ::FilePath::output_file:
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

string String::type_name() const {
    return "String";
}






} // closing namespace Type
} // closing namespace Input



