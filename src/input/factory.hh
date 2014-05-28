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


namespace Input {

using namespace std;


TYPEDEF_ERR_INFO( EI_KeyName, const string);
TYPEDEF_ERR_INFO( EI_TypeName, const string);
DECLARE_EXCEPTION( ExcNotRegistredClass, << "Key " << EI_KeyName::val
		<< " isn't registered in factory for type " << EI_TypeName::val << "!");

/*
template <class Type, class... Arguments>
class Creater {
public:
    shared_ptr<Type> operator()(Arguments ... args) {
        return make_shared<Type>(args...);
    }

};
*/


/**
 * This class implements more general factory mechanism to construct classes.
 *
 * One factory allows constructing derived classes of one base class. This class is determined
 * by template parameter Type and all descendants must implement constructor with same
 * parameters (given by template parameter Arguments). This constructor is called by factory.
 *
 * All descendants must contain:
 * 1. public static method for creation of new object stored to shared_ptr, this method is registered to factory
 * 2. constructor with parameters given by Arguments
 * 3. private static integer variable what is only for registration class to factory, this variable only allow
 *    to register class to factory and its implementation must call Factory::register_function what adds public
 *    static method to factory (see 1)
 *
 * Simple example of usage:
 @code
     class SomeDescendant : public SomeBase
     {
     public:
		/// create new object stored to shared pointer
		static std::shared_ptr< SomeBase > create_instance() {
			return std::make_shared< SomeDescendant >();
		}

		/// constructor
	    SomeDescendant() {}

     private:
     	/// registers class to factory
   	    static const int reg;
     }

     /// implementation of registration variable
     const int SomeDescendant::reg =
		 Input::Factory< SomeBase >::register_function("SomeDescendant", SomeDescendant::create_instance );
 @endcode
 *
 * Factory allow to accept constructor with one or more parameters. In this case Factory is also templated
 * by these parameters.
 * For example Factory< SomeBase, int, double > accepts constructors with two parameters (int, double).
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
   SomeBase * sb = a_rec.factory< SomeBase, double >(0.25);
 @endcode
 *
 *
 * TODO: used lambda function as second parameter of register_function method
 */
template <class Type, class... Arguments>
class Factory
{
public:
	/// Get the single instance of the factory
    static Factory * instance();

    /// Register lambda function that calls default constructor of Type.
    // Type of factory is BaseClass - not Type what is created
    //static int register_function(string class_name);

    /// register a factory function to create an instance of class_name
    static int register_function(string class_name, std::shared_ptr<Type>(* func)(Arguments...) );


    template <class Child>
    static int register_constructor(string class_name);


    /// create an instance of a registered class
    shared_ptr<Type> create(string name, Arguments... arguments);


private:
    /// a private constructor
    Factory(){}

    /// the registry of factory functions
    map<string, boost::any> factory_registry_;

};



} // closing namespace Input

// include implementation of templates and inline methods
#include "factory_impl.hh"

#endif // FIELD_RECORD_FACTORY_HH_
