// in the third_party/FADBAD++ dir, namespace "fadbad"

// undef asserts before including FADBAD
#ifdef ASSERT
    #undef ASSERT
#endif

#ifdef ASSERT_DBG
    #undef ASSERT_DBG
#endif

#include "fadbad.h"
#include "badiff.h"
#include "fadiff.h"
#include "tadiff.h"

// undef asserts from FADBAD
#ifdef ASSERT
    #undef ASSERT
#endif
#ifdef ASSERT_DBG
    #undef ASSERT_DBG
#endif

// define our asserts as in: 
// #include "system/asserts.hh"
#define ASSERT( expr) FEAL_ASSERT( expr)
#define ASSERT_DBG( expr) FEAL_ASSERT_DBG( expr)