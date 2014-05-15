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

/**
 * This class implements more general factory mechanism to construct classes given by TYPE key
 * of abstract record, without any explicit list of child classes.
 *
 * Factory can be used for every 'base class' that contains AbstractRecord and its descendants
 * are represented by Records derived from this AbstractRecord (e.g. FieldBase, ReactionBase etc.)
 *
 * Usage in base class:
 @code
	 class SomeBase
	 {
	 public:
		 SomeBase();

		 /// input abstract record
		 static Input::Type::AbstractRecord input_type;
	 };

	 Input::Type::AbstractRecord SomeBase::input_type
		 = Input::Type::AbstractRecord("SomeBase", "Abstract record desc.");


     /// Base class has to contain definition of factory
     template class Input::Factory< SomeBase >;
 @endcode
 *
 * Usage in derived class:
 @code
     class SomeDescendant : public SomeBase
     {
     public:
     	/// static variable for registration class to factory
   	    static const int reg;

   	    /// input record
   	    static Input::Type::Record input_type;

		/// static method for creation of new object stored to shared_ptr
		static std::shared_ptr< SomeBase > create_instance(int n_steps) {
			return std::make_shared< SomeDescendant >(n_steps);
		}

		/// constructor can have one ore more parameters
	    SomeDescendant(int n_steps) {}
     }

     /// input record is derived from abstract record of base class
     Input::Type::Record SomeDescendant::input_type =
     	 Input::Type::Record("SomeDescendant", "Record desc")
         .derive_from(SomeBase::input_type);

     /// implementation of registration variable
     const int SomeDescendant::reg =
		 Input::Factory< SomeBase >::register_function("SomeDescendant", SomeDescendant::create_instance );
 @endcode
 */
template <class Type>
class Factory
{
public:
    /// Get the single instance of the factory
    static Factory * instance();

    /// register a factory function to create an instance of class_name
    template <class... Arguments>
    static int register_function(string class_name, std::shared_ptr<Type>(* func)(Arguments...) );

    /// create an instance of a registered class
    template<class... Arguments>
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
