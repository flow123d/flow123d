/*
 * armadillo_setup.hh
 *
 *  Created on: Jul 4, 2016
 *      Author: jb
 */

#ifndef SRC_SYSTEM_ARMADILLO_SETUP_HH_
#define SRC_SYSTEM_ARMADILLO_SETUP_HH_


/**
 * This method sets particular ostream to armadillo which catch
 * errors reported by Armadillo and report Flow error message with full stacktrace.
 * This allows better determination of the source of the error.
 *
 * The call to the function may appear in the main() or in particular unit tests.
 */
void armadillo_setup();



#endif /* SRC_SYSTEM_ARMADILLO_SETUP_HH_ */
