/*
 * test_tools.hh
 *
 *  Created on: Jan 10, 2020
 *      Author: jb
 */

#ifndef TEST_TEST_TOOLS_HH_
#define TEST_TEST_TOOLS_HH_


#include "include/assert.hh"



#define ASSERT_THROW(expression, msg) {				            \
	bool success = false;										\
    try {                                                       \
        expression;                                             \
    } catch (const bparser::Exception & e) {                    \
    	std::string subwhat = std::string(e.what()).substr(0, std::string(msg).size()); \
		if (msg == subwhat ) \
			success = true;									    \
		else													\
			std::cout << "Wrong exception msg: " << e.what() << "\n";	\
    } catch (const std::exception &e) {                          \
    	throw e;												\
    }                                                           \
	if ( ! success) throw;								\
}


void fill_const(double *ptr, uint n, double v) {
	for(uint i=0; i< n; ++i) ptr[i] = v;
}

void fill_seq(double *ptr, uint n, double a, double st = 1) {
	for(uint i=0; i<n; a+=st, ++i) ptr[i] = a;
}

std::string print_vec(double *ptr, uint size) {
	std::stringstream s;
	s << "[";
	for(uint i=0; i<size; ++i)
		s << ptr[i] << ", " << "\n";
	s << "]";
	return s.str();
}

template< class T>
std::string print_vector(std::vector<T> x) {
	std::stringstream s;
	s << "[";
	for(uint i=0; i<x.size(); ++i)
		s << x[i] << ", " << "\n";
	s << "]";
	return s.str();
}

#endif /* TEST_TEST_TOOLS_HH_ */
