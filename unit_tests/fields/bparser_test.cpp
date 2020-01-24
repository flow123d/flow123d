/*
 * bparser_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include "include/parser.hh"
#include "./bparser_test_tools.hh"


TEST(BParser, all) {
	using namespace bparser;

	// Define own value vectors, preferably aligned.
	uint vec_size = 8;
	double v1[vec_size * 3];
	double v2[vec_size * 3];
	double vres[vec_size * 3];

	// Create parser, give the size of the value spaces.
	// That is maximal allocated space. Actual values and
	// active subset can be changed between evaluations of the expression.
	Parser p(vec_size);
	// parse an expression.
	p.parse("1 * v1 + cs1 * v2");

	// Get symbols used in the expressions.
	std::vector<std::string> symbols = p.symbols();
	std::cout << "Symbols: " << print_vector(p.symbols()) << "\n";
	// Set constants and variables

	// "cs1" constant with shape {}, i.e. scalar and values {2}.
	p.set_constant("cs1", {}, 	{2});
	// "cv1" vector constant with shape {3}
	p.set_constant("cv1", {3}, 	{1, 2, 3});
	// "v1" variable with shape {3}; v1 is pointer to the value space
	p.set_variable("v1", {3}, v1);
	p.set_variable("v2", {3}, v2);
	// Set the result variable (the last line of the expression)
	p.set_variable("_result_", {3}, vres);

	// Compile the expression into internall processor.
	p.compile();
	// Set arbitrary subset of the SIMD blocks in the maximal values space.
	// Here the full subset.
	p.set_subset({0, 1});
	// Evaluate
	p.run();
	// Result in the 'vres' value space.
	std::cout << print_vec(vres, 3*vec_size);
}
