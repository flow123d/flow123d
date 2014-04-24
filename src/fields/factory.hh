/*
 * field_record_factory.hh
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */

#ifndef FIELD_RECORD_FACTORY_HH_
#define FIELD_RECORD_FACTORY_HH_

#include "input/type_record.hh"

#include <memory>
#include <string>
#include <map>
#include <functional>

using namespace std;

// A helper class to register a factory function
class Registrar {
public:
    Registrar(string className, function<Input::Type::Record * (void)> classFactoryFunction);
};

// A preprocessor define used by derived classes
//#define REGISTER_CLASS(NAME, TYPE) static Registrar registrar(NAME, [](void) -> Input::Type::Record * { return TYPE::get_input_type();});

// The factory - implements singleton pattern!
class FieldRecordFactory
{
public:
    /// Get the single instance of the factory
    static FieldRecordFactory * instance();

    /// register a factory function to create an instance of className
    void register_factory_function(string name, function<Input::Type::Record*(void)> class_factory_function);

    /// create an instance of a registered class
    shared_ptr<Input::Type::Record> create(string name);

private:
    /// a private constructor
    FieldRecordFactory(){}

    /// the registry of factory functions
    map<string, function<Input::Type::Record*(void)>> field_factory_registry_;

};

#endif // FIELD_RECORD_FACTORY_HH_
