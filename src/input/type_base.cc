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
#include <boost/functional/hash.hpp>


#include "type_base.hh"
#include "type_record.hh"
#include "type_output.hh"
#include <boost/algorithm/string.hpp>


namespace Input {
namespace Type {

using namespace std;



/*******************************************************************
 * implementation of TypeBase
 */



TypeBase::TypeBase() {
    TypeBase::lazy_object_set().insert(this);
}



TypeBase::TypeBase(const TypeBase& other)
{
    TypeBase::lazy_object_set().insert(this);
}



TypeBase::~TypeBase() {
    TypeBase::LazyObjectsSet &set=TypeBase::lazy_object_set();
    TypeBase::LazyObjectsSet::iterator it =set.find(this);
    ASSERT( it != set.end(), "Missing pointer in lazy_object_set to '%s'.\n", this->type_name().c_str());
    TypeBase::lazy_object_set().erase(it);
}


bool TypeBase::is_valid_identifier(const string& key) {
  namespace ba = boost::algorithm;
  return ba::all( key, ba::is_lower() || ba::is_digit() || ba::is_any_of("_") );
}


string TypeBase::desc() const {
    stringstream ss;
    ss << OutputText(this,1);
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
    for (LazyTypeVector::iterator it=lazy_type_list().begin(); it!=lazy_type_list().end(); it++) {
        if (! (*it)->finish()) xprintf(PrgErr, "Can not finish '%s' during lazy_finish.\n", (*it)->type_name().c_str() );
    }

    lazy_type_list().clear();

}




TypeBase::LazyObjectsSet &TypeBase::lazy_object_set() {
    static LazyObjectsSet set_;
    return set_;
}



bool TypeBase::was_constructed(const TypeBase * ptr) {
    return lazy_object_set().find(ptr) != lazy_object_set().end();
}



std::ostream& operator<<(std::ostream& stream, const TypeBase& type) {
    return ( stream << OutputText(&type, 1) );
}



/**********************************************************************************
 * implementation of Type::Array
 */

std::size_t Array::content_hash() const
{
    std::size_t seed=0;
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, data_->lower_bound_);
    boost::hash_combine(seed, data_->upper_bound_);
    boost::hash_combine(seed, data_->type_of_values_->content_hash() );
    return seed;
}


bool Array::finish() {
	return data_->finish();
}



bool Array::ArrayData::finish()
{
	if (finished) return true;

	if (p_type_of_values != 0)
	{
	    if (! was_constructed(p_type_of_values) ) return false;

		if (dynamic_cast<const AbstractRecord *>(p_type_of_values) != 0)
		{
			AbstractRecord *ar = (AbstractRecord *)dynamic_cast<const AbstractRecord *>(p_type_of_values);
			boost::shared_ptr<const TypeBase> type_copy = boost::make_shared<const AbstractRecord>(*ar);
			type_of_values_ = type_copy;
			p_type_of_values = 0;
		}
		else if (dynamic_cast<const Record *>(p_type_of_values) != 0)
		{
			Record *r = (Record *)dynamic_cast<const Record *>(p_type_of_values);
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
		    xprintf(PrgErr, "Should not happen!\n");
	}

	return (finished = true);
}



string Array::type_name() const {
    return "array_of_" + data_->type_of_values_->type_name();
}

string Array::full_type_name() const {
	return "array_of_" + data_->type_of_values_->type_name();
}



bool Array::operator==(const TypeBase &other) const    {
    return  typeid(*this) == typeid(other) &&
              (*data_->type_of_values_ == static_cast<const Array *>(&other)->get_sub_type() );
}



bool Array::valid_default(const string &str) const {
    if ( this->match_size( 1 ) ) {
        return get_sub_type().valid_default( str );
    } else {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
    }
}


/**********************************************************************************
 * implementation and explicit instantiation of Array constructor template
 */

template <class ValueType>
Array::Array(const ValueType &type, unsigned int min_size, unsigned int max_size)
: data_(boost::make_shared<ArrayData>(min_size, max_size))
{
    // ASSERT MESSAGE: The type of declared keys has to be a class derived from TypeBase.
    BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, ValueType >::value) );
    ASSERT( min_size <= max_size, "Wrong limits for size of Input::Type::Array, min: %d, max: %d\n", min_size, max_size);

    // Records, AbstractRecords and Selections need not be initialized
    // at the moment, so we save the reference of type and update
    // the array later in finish().
    if ( (boost::is_base_of<Record, ValueType>::value ||
          boost::is_base_of<Selection, ValueType>::value)
         && ! TypeBase::was_constructed(&type) ) {
        data_->p_type_of_values = &type;
        TypeBase::lazy_type_list().push_back( boost::make_shared<Array>( *this ) );
    } else {
        data_->p_type_of_values = NULL;
        boost::shared_ptr<const TypeBase> type_copy = boost::make_shared<ValueType>(type);
        data_->type_of_values_ = type_copy;
        data_->finished=true;
    }
}

// explicit instantiation

#define ARRAY_CONSTRUCT(TYPE) \
template Array::Array(const TYPE &type, unsigned int min_size, unsigned int max_size)

ARRAY_CONSTRUCT(String);
ARRAY_CONSTRUCT(Integer);
ARRAY_CONSTRUCT(Double);
ARRAY_CONSTRUCT(Bool);
ARRAY_CONSTRUCT(FileName);
ARRAY_CONSTRUCT(Selection);
ARRAY_CONSTRUCT(Array);
ARRAY_CONSTRUCT(Record);
ARRAY_CONSTRUCT(AbstractRecord);


/**********************************************************************************
 * implementation of Type::Scalar ... and descendants.
 */

string Scalar::full_type_name() const {
    return type_name();
}

/**********************************************************************************
 * implementation of Type::Bool
 */


std::size_t Bool::content_hash() const
{
    std::size_t seed=0;
    boost::hash_combine(seed, type_name());
    return seed;
}


bool Bool::valid_default(const string &str) const {
    from_default(str);
    return true;
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


string Bool::type_name() const {
    return "Bool";
}


/**********************************************************************************
 * implementation of Type::Integer
 */

std::size_t Integer::content_hash() const
{
    std::size_t seed=0;
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, lower_bound_);
    boost::hash_combine(seed, upper_bound_);
    return seed;
}



bool Integer::match(std::int64_t value) const {
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



bool Integer::valid_default(const string &str) const
{
    from_default(str);
    return true;
}



string Integer::type_name() const {
    return "Integer";
}


/**********************************************************************************
 * implementation of Type::Double
 */


std::size_t Double::content_hash() const
{
    std::size_t seed=0;
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, lower_bound_);
    boost::hash_combine(seed, upper_bound_);
    return seed;
}



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



bool Double::valid_default(const string &str) const
{
    from_default(str);
    return true;
}




string Double::type_name() const {
    return "Double";
}


/**********************************************************************************
 * implementation of Type::FileName
 */

std::size_t FileName::content_hash() const
{
    std::size_t seed=0;
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, type_);
    return seed;
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


std::size_t String::content_hash() const
{
    std::size_t seed=0;
    boost::hash_combine(seed, type_name());
    return seed;
}



string String::type_name() const {
    return "String";
}




bool String::valid_default(const string &str) const {
    if (! match(str)) {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
    }
    return true;
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



