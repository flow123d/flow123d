/*
 * range_test.cpp
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest.hh>
#include <random>

#include "mesh/range_wrapper.hh"
#include "system/global_defs.h"


class TestHandler {
public:
    TestHandler(unsigned int size) {
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 20.99);

        values.resize(size);
        for (unsigned int i=0; i<size; ++i) values[i] = (unsigned int)dist(mt);
    }

    unsigned int value(unsigned int pos) const {
        return values[pos];
    }

private:
    std::vector<unsigned int> values;
};


class TestAccessor {
public:
    TestAccessor()
    : handler_(nullptr), pos_(0) {}

    TestAccessor(const TestHandler *handler, unsigned int pos)
    : handler_(handler), pos_(pos) {}

    unsigned int value() const {
    	return handler_->value(pos_);
    }

    bool operator==(const TestAccessor& other) {
    	return (pos_ == other.pos_);
    }

    operator unsigned int() const {
        return handler_->value(pos_);
    }

    inline void inc() {
        pos_++;
    }

private:
    const TestHandler *handler_;
    unsigned int pos_;

};


TEST(RangeWrapper, range) {
	const static unsigned int size = 20;
    TestHandler *test_handler = new TestHandler(size);
    std::vector<unsigned int> ret1, ret2;

    Range<TestAccessor> range(make_iter<TestAccessor>( TestAccessor(test_handler, 0) ),
                              make_iter<TestAccessor>( TestAccessor(test_handler, size) )
							 );
    for ( auto v : range ) ret1.push_back( v.value() );

    RangeConvert<TestAccessor, unsigned int> rangeC(make_iter<TestAccessor, unsigned int>( TestAccessor(test_handler, 0) ),
                                                    make_iter<TestAccessor, unsigned int>( TestAccessor(test_handler, size) )
							                       );
    for ( auto i : rangeC ) ret2.push_back( i );

    for (unsigned int i=0; i<size; ++i) {
        EXPECT_EQ( test_handler->value(i), ret1[i] );
        EXPECT_EQ( test_handler->value(i), ret2[i] );
    }

    delete test_handler;
}
