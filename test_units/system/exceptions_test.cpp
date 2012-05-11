/*
 * exceptions_test.cpp
 *
 *  Created on: May 10, 2012
 *      Author: jb
 */


#include "system/exceptions.hh"
#include <gtest/gtest.h>
#include <string>
#include <iostream>

/*
class SomeClass {

};

std::ostream & operator<<( std::ostream & stream, const SomeClass & sc)
{ return (stream << "Reporting SomeClass"); }
*/
//TYPEDEF_ERR_INFO( ExcStringValue, std::string );
TYPEDEF_ERR_INFO( ExcIntValue, int );
//TYPEDEF_ERR_INFO( ExcSomeClass, SomeClass );

DECLARE_EXCEPTION(ExcTest, << ExcIntValue_val <<"and :");
        //<< val(ExcStringValue::x) << "and a class:" << val<ExcSomeClass>() );




TEST(Exceptions, all) {
    //SomeClass sc;
    try {
        throw ExcTest() << ExcIntValue(42); // << ExcSomeClass(sc) << ExcStringValue("ahoj");
    } catch (ExceptionBase &e) {
        std::cout << e.what();
    }

}
