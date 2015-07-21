/*
 * type_generic.cc
 *
 *  Created on: May 1, 2012
 *      Author: jb
 */


#include <input/type_generic.hh>
#include <input/type_repository.hh>

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
    ASSERT(false, "Method valid_default can't be called for Parameter type.\n");
    return true;
}


TypeBase::MakeInstanceReturnType Parameter::make_instance(std::vector<ParameterPair> vec) const {
	ASSERT(false, "Method make_instance can't be called for type Parameter.\n");
	return std::make_pair( boost::make_shared<Parameter>(*this), this->parameter_map_ );
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
    ASSERT(false, "Method valid_default can't be called for Instance type.\n");
    return true;
}


const Instance &Instance::close() const {
	return *( Input::TypeRepository<Instance>::get_instance().add_type( *this ) );
}


// Implements @p TypeBase::make_instance.
TypeBase::MakeInstanceReturnType Instance::make_instance(std::vector<ParameterPair> vec) const {
	return generic_type_.make_instance(parameters_);
}

} // closing namespace Type
} // closing namespace Input
