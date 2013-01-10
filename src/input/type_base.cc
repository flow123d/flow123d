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

#include "system/system.hh"

#include <boost/type_traits.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>

#include "type_base.hh"
#include "type_record.hh"
#include <boost/algorithm/string.hpp>


namespace Input {
namespace Type {

using namespace std;




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



string TypeBase::desc() const {
    stringstream ss;
    reset_doc_flags();
    documentation(ss);
    return ss.str();
}



TypeBase::LazyTypeVector &TypeBase::lazy_type_list() {
    static LazyTypeVector lazy_type_list;
    return lazy_type_list;
}



void TypeBase::lazy_finish() {
    // TODO: dynamic cast as the switch may be expensive, in such case use some notification about type

    // first finish all lazy input types save Selection (we have to leave open Selection in AbstractType key TYPE)
    for (LazyTypeVector::iterator it=lazy_type_list().begin(); it!=lazy_type_list().end(); it++) {
        if (boost::dynamic_pointer_cast<Selection>(*it) == 0) {
            (*it)->finish();
        }
    }

    // then finalize abstract records so that no type can derive from them
    for (LazyTypeVector::iterator it=lazy_type_list().begin(); it!=lazy_type_list().end(); it++)
    {
        boost::shared_ptr<AbstractRecord> a_rec_ptr = boost::dynamic_pointer_cast<AbstractRecord>(*it);
        if ( a_rec_ptr!= 0) a_rec_ptr->no_more_descendants();
    }

    // at last finish all selections (including those in AbstractRecord)
    for (LazyTypeVector::iterator it=lazy_type_list().begin(); it!=lazy_type_list().end(); it++) (*it)->finish();

    lazy_type_list().clear();

}



std::ostream& operator<<(std::ostream& stream, const TypeBase& type) {
    type.reset_doc_flags();
    return type.documentation(stream);
}



/**********************************************************************************
 * implementation of Type::Array
 */

void Array::finish()
{
	empty_check();
	data_->finish();
}

void Array::ArrayData::finish()
{
	if (finished) return;

	if (p_type_of_values != 0)
	{
		if (dynamic_cast<const AbstractRecord *>(p_type_of_values) != 0)
		{
			AbstractRecord *ar = (AbstractRecord *)dynamic_cast<const AbstractRecord *>(p_type_of_values);
			ar->finish();
			boost::shared_ptr<const TypeBase> type_copy = boost::make_shared<const AbstractRecord>(*ar);
			type_of_values_ = type_copy;
			p_type_of_values = 0;
		}
		else if (dynamic_cast<const Record *>(p_type_of_values) != 0)
		{
			Record *r = (Record *)dynamic_cast<const Record *>(p_type_of_values);
			r->finish();
			boost::shared_ptr<const TypeBase> type_copy = boost::make_shared<const Record>(*r);
			type_of_values_ = type_copy;
			p_type_of_values = 0;
		}
		else if (dynamic_cast<const Selection *>(p_type_of_values) != 0)
		{
			Selection *s = (Selection *)dynamic_cast<const Selection *>(p_type_of_values);
			boost::shared_ptr<const TypeBase> type_copy = boost::make_shared<const Selection>(*s);
			type_of_values_ = type_copy;
			p_type_of_values = 0;
		}
		else if (dynamic_cast<const Array *>(p_type_of_values) != 0)
		{
			Array *a = (Array *)dynamic_cast<const Array *>(p_type_of_values);
			a->finish();
			boost::shared_ptr<const TypeBase> type_copy = boost::make_shared<const Array>(*a);
			type_of_values_ = type_copy;
			p_type_of_values = 0;
		}
	}

	finished = true;
}


std::ostream& Array::documentation(std::ostream& stream,DocType extensive, unsigned int pad) const {

	empty_check();

	switch (extensive) {
	case record_key:
		stream << "Array, size limits: [" << data_->lower_bound_ << ", " << data_->upper_bound_ << "] of type: " << endl;
		stream << setw(pad+4) << "";
		data_->type_of_values_->documentation(stream, record_key, pad+4);
		break;
	case full_after_record:
		data_->type_of_values_->documentation(stream, full_after_record, pad+4);
		break;
	case full_along:
		stream << "Array, size limits: [" << data_->lower_bound_ << ", " << data_->upper_bound_ << "] of type: " << endl;
		stream << setw(pad+4) << "";
		data_->type_of_values_->documentation(stream, record_key, pad+4);
		break;
	}

	return stream;
}



void  Array::reset_doc_flags() const {
	empty_check();
	data_->type_of_values_->reset_doc_flags();
}



string Array::type_name() const {
	empty_check();
    return "array_of_" + data_->type_of_values_->type_name();
}



bool Array::operator==(const TypeBase &other) const    {
	empty_check();
    return  typeid(*this) == typeid(other) &&
              (*data_->type_of_values_ == static_cast<const Array *>(&other)->get_sub_type() );
}



void Array::valid_default(const string &str) const {
    if ( this->match_size( 1 ) ) {
        get_sub_type().valid_default( str );
    } else {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
    }
}




/**********************************************************************************
 * implementation of Type::Scalar ... and descendants.
 */
void  Scalar::reset_doc_flags() const
{}

/**********************************************************************************
 * implementation of Type::Bool
 */

void Bool::valid_default(const string &str) const {
    from_default(str);
}



bool Bool::from_default(const string &str) const {
    if (str == "true" )  {
        return true;
    } else
    if (str == "false") {
        return false;
    } else {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
    }
}


std::ostream& Bool::documentation(std::ostream& stream,DocType extensive, unsigned int pad)  const {
    if (extensive == full_after_record) return stream;
    stream << "Bool";
    return stream;
}



string Bool::type_name() const {
    return "Bool";
}


/**********************************************************************************
 * implementation of Type::Integer
 */

bool Integer::match(int value) const {
    return ( value >=lower_bound_ && value <= upper_bound_);
}



int Integer::from_default(const string &str) const {
    std::istringstream stream(str);
    int value;
    stream >> value;

    if (stream && stream.eof() && match(value)) {
        return value;
    } else {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
    }
}



void Integer::valid_default(const string &str) const
{ from_default(str);}



std::ostream& Integer::documentation(std::ostream& stream,DocType extensive, unsigned int pad)  const {
    if (extensive == full_after_record) return stream;
    stream << "Integer in [" << lower_bound_ << ", " << upper_bound_ << "]";
    return stream;
}



string Integer::type_name() const {
    return "Integer";
}


/**********************************************************************************
 * implementation of Type::Double
 */

bool Double::match(double value) const {
    return ( value >=lower_bound_ && value <= upper_bound_);
}



double Double::from_default(const string &str) const {
    std::istringstream stream(str);
    double value;
    stream >> value;

    if (stream && stream.eof() && match(value)) {
        return value;
    } else {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
    }
}



void Double::valid_default(const string &str) const
{ from_default(str);}



std::ostream& Double::documentation(std::ostream& stream,DocType extensive, unsigned int pad)  const {
    if (extensive == full_after_record) return stream;
    stream << "Double in [" << lower_bound_ << ", " << upper_bound_ << "]";
    return stream;
}



string Double::type_name() const {
    return "Double";
}


/**********************************************************************************
 * implementation of Type::FileName
 */

std::ostream& FileName::documentation(std::ostream& stream,DocType extensive, unsigned int pad)  const {
    if (extensive == full_after_record) return stream;

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



bool FileName::match(const string &str) const {
    return (type_ == ::FilePath::input_file) || (str[0] != DIR_DELIMITER); // output files can not be absolute
}


/**********************************************************************************
 * implementation of Type::String
 */


std::ostream& String::documentation(std::ostream& stream,DocType extensive, unsigned int pad) const {

    if (extensive == full_after_record) return stream;
    stream << "String (generic)";
    return stream;
}



string String::type_name() const {
    return "String";
}




void String::valid_default(const string &str) const {
    if (! match(str)) {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
    }
}



string String::from_default(const string &str) const {
    valid_default(str);
    return str;
}



bool String::match(const string &str) const {
    return true;
}



} // closing namespace Type
} // closing namespace Input



