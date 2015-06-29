/*
 * type_generic.cc
 *
 *  Created on: May 1, 2012
 *      Author: jb
 */


#include <input/type_generic.hh>

#include <boost/functional/hash.hpp>



namespace Input {

namespace Type {


/*******************************************************************
 * implementation of Parameter
 */

Parameter::Parameter(const string & parameter_name)
: name_(parameter_name) {}


Parameter::Parameter(const Parameter & other)
: name_(other.name_) {}


string Parameter::type_name() const {
    return name_;
}


TypeBase::TypeHash Parameter::content_hash() const {
	TypeHash seed=0;
	boost::hash_combine(seed, "Parameter");
	boost::hash_combine(seed, type_name());

	return seed;
}


bool Parameter::valid_default(const string &str) const {
    if ( str != type_name() ) {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
    }
    return true;
}


/*******************************************************************
 * implementation of Instance
 */

Instance::Instance(const TypeBase &generic_type, std::vector<TypeBase::ParameterPair> parameters)
: generic_type_(generic_type), parameters_(parameters) {}


TypeBase::TypeHash Instance::content_hash() const {
	TypeHash seed=0;
	boost::hash_combine(seed, "Instance");
	boost::hash_combine(seed, generic_type_.content_hash() );
	for (std::vector<TypeBase::ParameterPair>::const_iterator it = parameters_.begin(); it!=parameters_.end(); it++) {
		boost::hash_combine(seed, (*it).first );
		boost::hash_combine(seed, (*it).second->content_hash() );
	}

	return seed;
}


bool Instance::valid_default(const string &str) const {
    return true;
}


bool Instance::finish() {
	// TODO returned type must be add to IST
	generic_type_.make_instance(parameters_);

	return true;
}


} // closing namespace Type
} // closing namespace Input
