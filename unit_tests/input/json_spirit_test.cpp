#include <flow_gtest.hh>

#include "input/json_spirit/json_spirit.h"


using namespace std;


/**
 * Tests if JSON string is valid, for invalid value returns
 * json_spirit::Error_position exception.
 *
 * JSON string must be wrapped in square brackets.
 */
void check_json_validity(std::string json_string) {
	json_spirit::mValue node;
	std::string wrapped_string("[ " + json_string + " ]");
	json_spirit::read_or_throw( wrapped_string, node );
}


TEST(JsonSpirit, parser) {
	// check JSON allowed inputs
	check_json_validity("true");
	check_json_validity("false");
	check_json_validity("10");
	check_json_validity("-3");
	check_json_validity("0.25");
	check_json_validity("-1.5");
	check_json_validity("\"ahoj\"");
	check_json_validity("\"some longer text\"");
	check_json_validity("[]");
	check_json_validity("[ 10, 20 ]");
	check_json_validity("[ \"val1\", \"val2\" ]");
	check_json_validity("{}");
	check_json_validity("{ key=10 }");
	check_json_validity("{ key1=5, key2=3.5, key3=\"ahoj\", key4=[ 0, 3, 5 ], key5={ a=1, b=\"str\" } }");

	// invalid JSON inputs
	EXPECT_THROW( { check_json_validity("ahoj"); }, json_spirit::Error_position);
	EXPECT_THROW( { check_json_validity("{ key=10, }"); }, json_spirit::Error_position);
	EXPECT_THROW( { check_json_validity("{ key= }"); }, json_spirit::Error_position);
	EXPECT_THROW( { check_json_validity("-3.6t5"); }, json_spirit::Error_position);
	EXPECT_THROW( { check_json_validity("1.23 4"); }, json_spirit::Error_position);
	EXPECT_THROW( { check_json_validity("1.23.4"); }, json_spirit::Error_position);
}
