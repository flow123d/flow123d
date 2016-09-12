/*
 * armadillo_setup.hh
 *
 *  Created on: Jul 4, 2016
 *      Author: jb
 */

#ifndef SRC_SYSTEM_ARMADILLO_TOOLS_HH_
#define SRC_SYSTEM_ARMADILLO_TOOLS_HH_

#include <string>
#include <armadillo>

/**
 * This method sets particular ostream to armadillo which catch
 * errors reported by Armadillo and report Flow error message with full stacktrace.
 * This allows better determination of the source of the error.
 *
 * The call to the function may appear in the main() or in particular unit tests.
 */
void armadillo_setup();


/**
 * Format field_value (i.e. matrix, vector scalar of double or int) into YAML string.
 */
template<class T>
std::string field_value_to_yaml(const T &mat, unsigned int prec = 5);


#endif /* SRC_SYSTEM_ARMADILLO_TOOLS_HH_ */
