/*
 * check_error_test.cpp
 *
 *  Created on: May 24, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>
#include "system/exceptions.hh"
#include "system/global_defs.h"


TEST(CheckError, error_message) {
	unsigned int err_code;

	err_code = 0;
	CHKERR( err_code );

	err_code = 1;
	EXPECT_THROW_WHAT( { CHKERR( err_code ); }, ExcChkErr, "1" );
}

TEST(CheckError, assert_message) {
	unsigned int err_code;

	err_code = 0;
	CHKERR_ASSERT( err_code );

	err_code = 1;
	EXPECT_THROW_WHAT( { CHKERR_ASSERT( err_code ); }, ExcChkErrAssert, "1" );
}
