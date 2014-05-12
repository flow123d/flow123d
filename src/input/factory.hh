/*
 * field_record_factory.hh
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */

#ifndef FIELD_RECORD_FACTORY_HH_
#define FIELD_RECORD_FACTORY_HH_

#include "input/accessors.hh"
#include "input/type_record.hh"
//#include "fields/field_base.hh"

#include <memory>
#include <string>
#include <map>
#include <functional>
#include <boost/functional/factory.hpp>
#include <boost/any.hpp>


namespace Input {

using namespace std;

// The factory - implements singleton pattern!
template <class Type>
class Factory
{
	friend class AbstractRecord;
public:
    /// Get the single instance of the factory
    static Factory * instance();

    /// register a factory function to create an instance of className
    void register_function(string name, const boost::any& func);

    static void add(string name, const boost::any& func) {};

    /// create an instance of a registered class
    template<class... Arguments>
    shared_ptr<Type> create(string name, Arguments... arguments);


private:
    /// a private constructor
    Factory(){}

    /// the registry of factory functions
    map<string, boost::any> factory_registry_;

};


// A helper class to register a factory function
template<class Type>
class Registrar {
public:
	typedef typename Type::FactoryBaseType BaseType;

	Registrar(string className);

	//Registrar(string className, const boost::any& func);
	template <class... args>
	Registrar(string class_name, std::shared_ptr<BaseType>(* func)(args...) ) {
	    auto func_wrapper = std::function<std::shared_ptr<BaseType>(args...)>(func);
	    Factory<BaseType>::instance()->register_function(class_name, func_wrapper);
	}
};


} // closing namespace Input

// include implementation of templates and inline methods
#include "factory_impl.hh"

#endif // FIELD_RECORD_FACTORY_HH_
