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
    ASSERT(false, "Method valid_default can't be called for Parameter '%s'.\n", this->name_.c_str());
    return true;
}


TypeBase::MakeInstanceReturnType Parameter::make_instance(std::vector<ParameterPair> vec) const {
	ParameterMap parameter_map;
	for (std::vector<ParameterPair>::iterator vec_it=vec.begin(); vec_it!=vec.end(); vec_it++) {
		if ( (*vec_it).first == this->name_ ) {
			parameter_map[(*vec_it).first] = (*vec_it).second->content_hash();
			return std::make_pair( (*vec_it).second, parameter_map );
		}
	}
    THROW( ExcParamaterNotSubsituted() << EI_Object(this->name_));
	return std::make_pair( boost::make_shared<Parameter>(*this), parameter_map );
}


bool Parameter::finish(bool is_generic) {
	if (!is_generic) THROW( ExcParamaterInIst() << EI_Object(this->name_));
	return true;
}


/*******************************************************************
 * implementation of Instance
 */

Instance::Instance(TypeBase &generic_type, std::vector<TypeBase::ParameterPair> parameters)
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


bool Instance::finish(bool is_generic) {
	return generic_type_.finish(true);
}


std::string print_parameter_vec(std::vector<TypeBase::ParameterPair> vec) {
	stringstream ss;
	for (std::vector<TypeBase::ParameterPair>::const_iterator vec_it = vec.begin(); vec_it!=vec.end(); vec_it++) {
		if (vec_it != vec.begin()) ss << ", ";
		ss << "'" << vec_it->first << "'" << endl;
	}

	return ss.str();
}


// Implements @p TypeBase::make_instance.
TypeBase::MakeInstanceReturnType Instance::make_instance(std::vector<ParameterPair> vec) const {
	// check if instance is created
	if (created_instance_.first) {
		return created_instance_;
	}

	try {
		created_instance_ = generic_type_.make_instance(parameters_);
	} catch (ExcParamaterNotSubsituted &e) {
        e << EI_ParameterList( print_parameter_vec(parameters_) );
        throw;
	}
#ifdef FLOW123D_DEBUG_ASSERTS
	for (std::vector<TypeBase::ParameterPair>::const_iterator vec_it = parameters_.begin(); vec_it!=parameters_.end(); vec_it++) {
		ParameterMap::iterator map_it = created_instance_.second.find( vec_it->first );
		ASSERT(map_it != created_instance_.second.end(), "Unused parameter '%s' in input type instance with parameters: %s.\n",
				vec_it->first.c_str(), print_parameter_vec(parameters_).c_str());
	}
#endif
	return created_instance_;
}

} // closing namespace Type
} // closing namespace Input
