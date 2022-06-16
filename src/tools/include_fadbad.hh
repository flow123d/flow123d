// in the third_party/FADBAD++ dir, namespace "fadbad"

// undef asserts before including FADBAD
#ifdef ASSERT
    #undef ASSERT
#endif

#ifdef ASSERT_PERMANENT
    #undef ASSERT_PERMANENT
#endif

#include "fadbad.h"
#include "badiff.h"
#include "fadiff.h"
#include "tadiff.h"

// undef asserts from FADBAD
#ifdef ASSERT
    #undef ASSERT
#endif
#ifdef ASSERT_PERMANENT
    #undef ASSERT_PERMANENT
#endif

// define our asserts as in: 
// #include "system/asserts.hh"
#define ASSERT( expr) FEAL_ASSERT( expr)
#define ASSERT_PERMANENT( expr) FEAL_ASSERT_PERMANENT( expr)
