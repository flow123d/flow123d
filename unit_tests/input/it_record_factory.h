/*
 * it_record_factory.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef IT_RECORD_FACTORY_HH_
#define IT_RECORD_FACTORY_HH_

#include <input/type_record.hh>


class ITRecordRegistrar
{
public:
	ITRecordRegistrar(string name, string desc);
};


class ITAbstractRecordRegistrar
{
public:
	ITAbstractRecordRegistrar(string name, string desc);
};


class ITRecordFactory
{
public:
	static ITRecordFactory * instance()
	{
		static ITRecordFactory factory_;
		return &factory_;
	}

	std::shared_ptr<Input::Type::Record> get_record(string name)
	{
		Input::Type::Record * instance = nullptr;

		// find name in the registry and call factory method.
		auto it = factoryRecordRegistry.find(name);
		if(it != factoryRecordRegistry.end())
			instance = it->second;

		// wrap instance in a shared ptr and return
		if(instance != nullptr)
			return std::shared_ptr<Input::Type::Record>(instance);
		else
			return nullptr;
	}

	std::shared_ptr<Input::Type::AbstractRecord> get_abstract_record(string name)
	{
		Input::Type::AbstractRecord * instance = nullptr;

		// find name in the registry and call factory method.
		auto it = factoryAbstractRecordRegistry.find(name);
		if(it != factoryAbstractRecordRegistry.end())
			instance = it->second;

		// wrap instance in a shared ptr and return
		if(instance != nullptr)
			return std::shared_ptr<Input::Type::AbstractRecord>(instance);
		else
			return nullptr;
	}

	void register_factory_record(string name, Input::Type::Record * rec)
	{
		factoryRecordRegistry[name] = rec;
	}

	void register_factory_abstract_record(string name, Input::Type::AbstractRecord * rec)
	{
		factoryAbstractRecordRegistry[name] = rec;
	}

private:
	/// default constructor
	ITRecordFactory() {};

	std::map<string, Input::Type::Record *> factoryRecordRegistry;
	std::map<string, Input::Type::AbstractRecord *> factoryAbstractRecordRegistry;
};

ITRecordRegistrar::ITRecordRegistrar(string name, string desc)
{
	//cout << "ITRecordRegistrar constructor: " << name << endl;
	ITRecordFactory::instance()->register_factory_record(name, new Input::Type::Record(name, desc));
}

ITAbstractRecordRegistrar::ITAbstractRecordRegistrar(string name, string desc)
{
	//cout << "ITAbstractRecordRegistrar constructor: " << name << endl;
	ITRecordFactory::instance()->register_factory_abstract_record(name, new Input::Type::AbstractRecord(name, desc));
}

#endif /* IT_RECORD_FACTORY_HH_ */
