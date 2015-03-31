/*
 * field_record_factory.hh
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */

#ifndef FIELD_RECORD_FACTORY_HH_
#define FIELD_RECORD_FACTORY_HH_

#include <memory>
#include <string>
#include <map>
#include <functional>
#include <boost/functional/factory.hpp>
#include <boost/any.hpp>

#include "system/exceptions.hh"


namespace Input {

using namespace std;


TYPEDEF_ERR_INFO( EI_KeyName, const string);
TYPEDEF_ERR_INFO( EI_TypeName, const string);
DECLARE_EXCEPTION( ExcNotRegistredClass, << "Key " << EI_KeyName::val
		<< " isn't registered in factory for type " << EI_TypeName::val << "!");


/**
 * This class implements more general factory mechanism to construct classes.
 *
 * One factory allows constructing derived classes of one base class. This class is determined
 * by template parameter Type and all descendants must implement constructor with same
 * parameters (given by template parameter Arguments). This constructor is called by factory.
 *
 * All descendants must contain:
 * 1. constructor with parameters given by Arguments
 * 2. declaration of parent class as typedef with name FactoryBaseType
 * 3. private static integer variable what is only for registration class to factory, this variable only allow
 *    to register class to factory and its implementation must call Factory::register_class what adds constructor
 *    of class to factory
 *
 * Simple example of usage:
 @code
     class SomeDescendant : public SomeBase
     {
     public:
        /// typedef of parent class
        typedef SomeBase FactoryBaseType;

		/// constructor
	    SomeDescendant() {}

     private:
     	/// registers class to factory
   	    static const int reg;
     }

     /// implementation of registration variable
     const int SomeDescendant::reg =
		 Input::register_class< SomeDescendant >("SomeDescendant");
 @endcode
 *
 * Factory allow to accept constructor with one or more parameters. In this case Factory is
 * also templated by these parameters. For example Factory< SomeBase, int, double > accepts
 * constructors with two parameters (int, double).
 *
 * If registered class is templated the following design have to be used:
 @code
     /// Example of class templated by integer parameter
     template <int dimension>
     class SomeDescendant : public SomeBase<dimension>
     {
     public:

		/// constructor
	    SomeDescendant(double time) {}
  	    ...
     }

     /// implementation of registration variable
     const int SomeDescendant::reg =
		 Input::register_class< SomeDescendant<dimension>, double >("SomeDescendant");
 @endcode
 *
 * Factory can be used in two ways:
 * - through Factory::create method
 *   Example for constructor with one parameter:
 @code
     SomeBase * sb = Input::Factory< SomeBase, double >::instance()->create("SomeDescendant", 0.1);
 @endcode
 * - through AbstractRecord::factory method. This possibility can be used if base class has defined
 *   AbstractRecord and its descendants contain Record derived from this AbstractRecord.
 *   Example for same constructor:
 @code
     AbstractRecord a_rec = record.val<AbstractRecord>("problem");
     SomeBase * sb = a_rec.factory< SomeBase, double>(0.25);
     // If arguments types can be deduced (by compiler) from actual arguments one can even use:
     SomeBase * sb = a_rec.factory< SomeBase >(0.25);
 @endcode
 *
 */
template <class Type, class... Arguments>
class Factory
{
public:
	/// Get the single instance of the factory
    static Factory * instance();


    /// Register lambda function that calls default constructor of Type.
    template <class Child>
    static int register_class(string class_name);


    /// create an instance of a registered class
    shared_ptr<Type> const create(string name, Arguments... arguments) const;


private:
    /// a private constructor
    Factory(){}

    /// the registry of factory functions
    map<string, boost::any> factory_registry_;

};


/**
 * Function allows simplified call of registering class to factory.
 *
 * It is used for declaration of registration variable.
 * @see Factory
 *
 * Example of usage:
 @code
     const int SomeClass::reg =
		 Input::register_class< SomeClass >("SomeClass");
 @endcode
 */
template <class ChildType, class... Arguments>
int register_class(string class_name);


} // closing namespace Input

// include implementation of templates and inline methods
#include "factory_impl.hh"

#endif // FIELD_RECORD_FACTORY_HH_
