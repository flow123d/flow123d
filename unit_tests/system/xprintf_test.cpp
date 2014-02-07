/*
 * xprintf_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>
#include "system/system.hh"
#include "system/exceptions.hh"


TEST(Xprintf, error_messages) {
	std::string key = "blue";
	std::string type = "ColorSelection";

	EXPECT_THROW_WHAT(
			{ xprintf(Err, "Unknown solver type. Internal error.\n"); },
			ExcXprintfMsg,
			"Unknown solver type. Internal error."
			);

	EXPECT_THROW_WHAT(
			{ xprintf(PrgErr, "Name '%s' already exists in Selection: %s\n", key.c_str(), type.c_str()); },
			ExcXprintfMsg,
			"already exists in Selection:"
			);

	EXPECT_THROW_WHAT(
			{ xprintf(UsrErr, "Unsupported problem type: %d\n", 10); },
			ExcXprintfMsg,
			"Unsupported problem type:"
			);
}
