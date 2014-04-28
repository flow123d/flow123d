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
#include "fields/field_base.hh"

#include <memory>
#include <string>
#include <map>
#include <functional>

namespace Input {

using namespace std;

// A helper class to register a factory function
template<class Type>
class Registrar {
public:
    Registrar(string className);
};


// The factory - implements singleton pattern!
class Factory
{
	friend class AbstractRecord;
public:
    /// Get the single instance of the factory
    static Factory * instance();

//protected:
    /// register a factory function to create an instance of className
    template<class Type>
    void register_function(string name, function<Type*(void)> class_factory_function);

    /// create an instance of a registered class
    template<class Type>
    shared_ptr<Type> create(AbstractRecord &rec);

private:
    /// a private constructor
    Factory(){}

    /// the registry of factory functions
    map<string, function<void*(void)>> field_factory_registry_;

};


} // closing namespace Input

// include implementation of templates and inline methods
#include "factory_impl.hh"

#endif // FIELD_RECORD_FACTORY_HH_
