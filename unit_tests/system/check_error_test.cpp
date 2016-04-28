/*
 * check_error_test.cpp
 *
 *  Created on: May 24, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>
#include "system/exceptions.hh"
#include "system/global_defs.h"
#include "system/system.hh"


TEST(CheckError, error_message) {
	unsigned int err_code;

	err_code = 0;
	chkerr( err_code );

	err_code = 1;
	EXPECT_THROW_WHAT( { chkerr( err_code ); }, ExcChkErr, "1" );
}

TEST(CheckError, assert_message) {
	unsigned int err_code;

	err_code = 0;
	chkerr_assert( err_code );

#ifdef FLOW123D_DEBUG_ASSERTS
    err_code = 1;
    EXPECT_THROW_WHAT( { chkerr_assert( err_code ); }, ExcChkErr, "1" );
#endif
}

#ifdef FLOW123D_DEBUG_ASSERTS
TEST(ASSERTS, assert_ptr) {
    void * test_ptr = nullptr;
    EXPECT_THROW_WHAT( {ASSERT_PTR(test_ptr);}, ExcAssertMsg, "test_ptr");
}
#endif
