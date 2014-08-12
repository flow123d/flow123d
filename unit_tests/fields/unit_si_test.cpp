/*
 * units_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>

#include "fields/unit.hh"


TEST(Units, printout) {
	Unit unit_1 = Unit().m(2).kg(1).s(-2);
	std::cout << unit_1.print() << std::endl;
	EXPECT_STREQ(std::string("m^2.kg.s^(-2)").c_str(), unit_1.print().c_str());

	Unit unit_2 = Unit().m(1).kg(0).s(-2);
	std::cout << unit_2.print() << std::endl;
	EXPECT_STREQ(std::string("m.s^(-2)").c_str(), unit_2.print().c_str());
}
