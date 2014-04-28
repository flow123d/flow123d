/*
 * field_record_factory.cc
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */


#include "fields/field_base.hh"
#include "input/factory.hh"


namespace Input {

using namespace std;


Factory * Factory::instance()
{
    static Factory factory;
    return &factory;
}

} // closing namespace Input
