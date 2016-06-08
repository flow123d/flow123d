/*
 * exceptions_test.cpp
 *
 *  Created on: May 10, 2012
 *      Author: jb
 */

#define FLOW123D_DEBUG_ASSERTS_WITHOUT_MPI
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include "system/exceptions.hh"

#include <string>
#include <iostream>
#include <typeinfo>


// Using a class with overloaded "<<" operator inside exception's messages.
class SomeClass {
};
std::ostream & operator<<( std::ostream & stream, const SomeClass & sc)
{ return (stream << "Reporting SomeClass"); }



// Declare ERR_INFO types and exceptions to test.

TYPEDEF_ERR_INFO( EI_StringValue, std::string );
TYPEDEF_ERR_INFO( EI_IntValue, int );
TYPEDEF_ERR_INFO( EI_SomeClass, SomeClass );

DECLARE_EXCEPTION(ExcInt, << "Int exception value: "<< EI_IntValue::val );
DECLARE_EXCEPTION(ExcString, << "Str exception :" << EI_StringValue::qval );
DECLARE_EXCEPTION(ExcClass, << "Class exception :" << EI_SomeClass::qval );
DECLARE_FATAL_EXCEPTION(ExcFatal, << "Test exception :" << EI_StringValue::qval );



// Test exception with EI_::qval - quoted value
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

// Test exception with EI_ without given value
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

// Test exception with EI_::val
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


// Test EI_ for passing a general class.
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


// Test EI_ for passing a fatal exception.
TEST(Exceptions, Fatal) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    try {
        THROW( ExcFatal() << EI_StringValue("ahoj") );
    } catch (ExceptionBase &e) {
        std::ostringstream ss;
        e.print_info(ss);
        EXPECT_EQ(std::string("Test exception :'ahoj'\n"), ss.str());

        std::cerr << e.what();
    }
}




//------------------------------------------------------------------------
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
//----------------------------------------------------------------------------

// Test exception and ErrorInfo declared inside class inside namespace.
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



//------------------------------------------------------------------------
class BlackBox {
public:
    inline int eval() const {return 42;}
};

TYPEDEF_ERR_INFO( EI_BlackBox, BlackBox );
DECLARE_EXCEPTION(ExcStorageTypeMismatch, << "stored is value of type " << EI_BlackBox::ref(*this).eval() );
//----------------------------------------------------------------------------

// Test calling a method to get error message from ErrorInfo object. Using EI_ref(_exc).
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




//------------------------------------------------------------------------
// Test stack trace
void zero()  { THROW(ExcString() << EI_StringValue("BOOM!")); }
void one()   { zero(); }
void two()   { one(); }
void three() { two(); }
TEST(Exceptions, stack_trace) {
    try {
        three();
    } catch (ExcString &e) {
        std::cerr << e.what();
    }

}

#include "system/global_defs.h"
//------------------------------------------------------------------------
// Test ASSERT using exceptions.
TEST(Exceptions, assert_msg) {

    int zero=1;
#ifdef FLOW123D_DEBUG_ASSERTS
    EXPECT_THROW_WHAT( {OLD_ASSERT(zero==0, "Zero is %d not zero.\n", zero);} , ExcAssertMsg, "Violated Assert! Zero is 1 not zero." );
#endif

}


// Empty class. used for ASSERT_PTR test
class EmptyObj {};


// Test of new asserts.
TEST(FealAssert, assert) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    std::string s1 = "feal";
    std::string s2 = "assert";

    try {
        FEAL_ASSERT(s1.empty() && s2.empty())(s1)(s2).error();
    } catch (feal::Exc_assert &e) {
        std::cout << e.what();
    }

    EXPECT_DEATH( {FEAL_ASSERT(s1.empty() && s2.empty())(s1)(s2);}, ".*");

    // only in debug mode
    try {
        FEAL_ASSERT_DBG(s1.empty() && s2.empty())(s1)(s2).error();
    } catch (feal::Exc_assert &e) {
        std::cout << e.what();
    }

    // comparative asserts
    {
    	int i=5, j=4;
    	EmptyObj *empty = nullptr;

    	EXPECT_THROW_WHAT( { ASSERT_LT(i, j).error(); }, feal::Exc_assert, "Expression: 'i < j'" );
    	EXPECT_THROW_WHAT( { ASSERT_LE(i, j).error(); }, feal::Exc_assert, "Expression: 'i <= j'" );
    	EXPECT_THROW_WHAT( { ASSERT_GT(j, i).error(); }, feal::Exc_assert, "Expression: 'j > i'" );
    	EXPECT_THROW_WHAT( { ASSERT_GE(j, i).error(); }, feal::Exc_assert, "Expression: 'j >= i'" );
    	EXPECT_THROW_WHAT( { ASSERT_EQ(i, j).error(); }, feal::Exc_assert, "Expression: 'i == j'" );
    	EXPECT_THROW( { ASSERT_PTR(empty).error(); }, feal::Exc_assert );
    }

}
