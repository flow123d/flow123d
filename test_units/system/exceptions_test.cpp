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
#include <typeinfo>

#include <boost/lambda/lambda.hpp>

class SomeClass {

};

std::ostream & operator<<( std::ostream & stream, const SomeClass & sc)
{ return (stream << "Reporting SomeClass"); }

TYPEDEF_ERR_INFO( EI_StringValue, std::string );
TYPEDEF_ERR_INFO( EI_IntValue, int );
TYPEDEF_ERR_INFO( EI_SomeClass, SomeClass );

DECLARE_EXCEPTION(ExcInt, << "Int exception value: "<< EI_IntValue::val );
DECLARE_EXCEPTION(ExcString, << "Str exception :" << EI_StringValue::qval );
DECLARE_EXCEPTION(ExcClass, << "Class exception :" << EI_SomeClass::qval );



TEST(Exceptions, String) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    try {
        THROW( ExcString() << EI_StringValue("ahoj") );
    } catch (ExceptionBase &e) {
        std::ostringstream ss;
        e.print_info(ss);
        EXPECT_EQ(std::string("Str exception :'ahoj'\n"), ss.str());

        std::cerr << e.what();
    }
}

TEST(Exceptions, no_value) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    try {
        THROW( ExcString() );
    } catch (ExceptionBase &e) {
        std::ostringstream ss;
        e.print_info(ss);
        EXPECT_EQ(std::string("Str exception :NO_VALUE\n"), ss.str());

        std::cerr << e.what();
    }
}

TEST(Exceptions, Int) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    try {
        THROW( ExcInt() << EI_IntValue(42) );
    } catch (ExceptionBase &e) {
        std::ostringstream ss;
        e.print_info(ss);
        EXPECT_EQ(std::string("Int exception value: 42\n"), ss.str());

        std::cerr << e.what();
    }
}



TEST(Exceptions, Class) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    try {
        SomeClass sc;
        THROW( ExcClass() << EI_SomeClass(sc) );
    } catch (ExceptionBase &e) {
        std::ostringstream ss;
        e.print_info(ss);
        EXPECT_EQ(std::string("Class exception :'Reporting SomeClass'\n"), ss.str());

        std::cout << e.what();
    }
}





namespace Inner {
class SomeClass {
public:
    TYPEDEF_ERR_INFO(EI_SomeClass, SomeClass);
    DECLARE_EXCEPTION(ExcSomeClass, << EI_SomeClass::val);

    void test() {
        THROW( ExcSomeClass() << EI_SomeClass(*this) );
    }
};

std::ostream & operator<<( std::ostream & stream, const SomeClass & sc)
{ return (stream << "Reporting Inner::SomeClass"); }

} // namespace Inner

TEST(Exceptions, InnerClass) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    try {
        Inner::SomeClass sc;
        sc.test();

    } catch (ExceptionBase &e) {
        std::ostringstream ss;
        e.print_info(ss);
        EXPECT_EQ(std::string("Reporting Inner::SomeClass\n"), ss.str());

        std::cout << e.what();
    }
}

class BlackBox {
public:
    inline int eval() const {return 42;}
};

TYPEDEF_ERR_INFO( EI_BlackBox, BlackBox );
DECLARE_EXCEPTION(ExcStorageTypeMismatch, << "stored is value of type " << EI_BlackBox::ref(_exc).eval() );


TEST(Exceptions, GettingValue) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    try {
        BlackBox bb;
        THROW( ExcStorageTypeMismatch() << EI_BlackBox(bb) );
    } catch (ExcStorageTypeMismatch &e) {
        std::ostringstream ss;
        e.print_info(ss);
        EXPECT_EQ(std::string("stored is value of type 42\n"), ss.str());

        std::cerr << e.what();
    }

}
