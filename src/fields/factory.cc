/*
 * field_record_factory.cc
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */


#include "field_record_factory.hh"


using namespace std;


Registrar::Registrar(string name, function<Input::Type::Record * (void)> class_factory_function)
{
    // register the class factory function
	FieldRecordFactory::instance()->register_factory_function(name, class_factory_function);
}


FieldRecordFactory * FieldRecordFactory::instance()
{
    static FieldRecordFactory factory;
    return &factory;
}


void FieldRecordFactory::register_factory_function(string name, function<Input::Type::Record*(void)> class_factory_function)
{
    // register the class factory function
	field_factory_registry_[name] = class_factory_function;
}


shared_ptr<Input::Type::Record> FieldRecordFactory::create(string name)
{
	Input::Type::Record * instance = nullptr;

    // find name in the registry and call factory method.
    auto it = field_factory_registry_.find(name);
    if(it != field_factory_registry_.end())
        instance = it->second();

    // wrap instance in a shared ptr and return
    if(instance != nullptr)
        return std::shared_ptr<Input::Type::Record>(instance);
    else
        return nullptr;
}
