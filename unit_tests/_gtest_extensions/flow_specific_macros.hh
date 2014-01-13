/*
 * flow_specific_macros.hh
 *
 *  Created on: Dec 16, 2013
 *      Author: jb
 */

#ifndef FLOW_SPECIFIC_MACROS_HH_
#define FLOW_SPECIFIC_MACROS_HH_

/**
 * Macro to test ASSERTS. It is similar to EXPECT_THROW_WHAT but
 * - it sets exception type to ExcAssertMsg
 * - it is defined empty if ASSERTS are off
 */
#ifdef DEBUG_ASSERTS
#define	EXPECT_ASSERT_DEATH(statement, pattern) EXPECT_THROW_WHAT(statement,ExcAssertMsg, pattern)
#else
#define	EXPECT_ASSERT_DEATH(statement, pattern)
#endif

/**
 * Possibility to turn off all death tests. Eg. when running under valgrind.
 * ... they are fragile and usually cause SEGFAULT
 */
#ifdef UNIT_TESTS_NO_DEATH_TESTS

#undef EXPECT_ASSERT_DEATH
#define EXPECT_ASSERT_DEATH(statement, pattern)

#undef EXPECT_THROW_WHAT
#define EXPECT_THROW_WHAT(statement, exception, pattern)

#undef EXPECT_THROW
#define EXPECT_THROW(statement, exception)

#undef EXPECT_DEATH
#define EXPECT_DEATH(statement, pattern)

#undef ASSERT_THROW
#define ASSERT_THROW(statement, exception)

#undef ASSERT_DEATH
#define ASSERT_DEATH(statement, pattern)
#endif




#endif /* FLOW_SPECIFIC_MACROS_HH_ */
