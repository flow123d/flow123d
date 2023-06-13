/*
 * bparser_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <string>
#include "assert.hh"
#include "parser.hh"


TEST(BParser, all) {
    using namespace bparser;

	// Define own value vectors, preferably aligned.
	uint vec_size = 24;
	double coords[vec_size * 3];
	fill_seq(coords, 3 * vec_size, 0.25, 0.5);
	double *x = coords;
	double *y = coords + vec_size;
	double *z = coords + 2*vec_size;
	double vres[vec_size];
	fill_const(vres, vec_size, -1);

	// Create parser, give the size of the value spaces.
	// That is maximal allocated space. Actual values and
	// active subset can be changed between evaluations of the expression.
	Parser p(vec_size);
	p.parse("1 * x + cs1 * y + x * z");

	std::cout << "Variables: " << print_vector(p.free_symbols()) << "\n";

	// Set constants and variables
	// "cs1" constant with shape {}, i.e. scalar and values {2}.
	p.set_constant("cs1", {}, 	{2});
	// "cv1" vector constant with shape {3}
	p.set_constant("cv1", {3}, 	{1, 2, 3});
	// "v1" variable with shape {3}; v1 is pointer to the value space
	p.set_variable("x", {}, x);
	p.set_variable("y", {}, y);
	p.set_variable("z", {}, z);
	// Set the result variable (the last line of the expression)
	p.set_variable("_result_", {}, vres);

	std::cout << "Symbols: " << print_vector(p.symbols()) << "\n";

	// Compile the expression into internall processor.
	p.compile();
	// Set arbitrary subset of the SIMD blocks in the maximal values space.
	// Here the full subset.
	uint simd_size = bparser::get_simd_size();
	std::vector<uint> ss = std::vector<uint>(vec_size/simd_size);
	for (uint i = 0; i < vec_size/simd_size; i++){
		ss[i] = i;
	}
	p.set_subset(ss);
	// Evaluate
	p.run();
	// Result in the 'vres' value space.
	std::cout << print_vec(vres, vec_size); //*/
}
