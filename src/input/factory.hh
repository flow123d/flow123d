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

// The factory - implements singleton pattern!
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
