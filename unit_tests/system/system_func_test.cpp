/*
 * system_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */



#include <gtest/gtest.h>
#include <sstream>
#include <string>
#include "system/system.hh"


string input = R"CODE(

section A
AAA
  // section B
  section B
BBB


)CODE";

TEST(SystemFunctions, skip_to) {
    std::stringstream ss(input);

    string line;
    skip_to(ss, "section A");
    ss >> line;
    EXPECT_EQ("AAA", line);
    skip_to(ss, "section B");
    std::getline(ss, line);
    EXPECT_EQ("BBB", line);

}
