/*
 * unit_si_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>

#include "fields/unit_si.hh"


TEST(UnitSI, printout) {
	UnitSI unit_1 = UnitSI().m(2).kg(1).s(-2);
	std::cout << unit_1.print() << std::endl;
	EXPECT_STREQ(std::string("$[m^{2}kgs^{-2}]$").c_str(), unit_1.print().c_str());

	UnitSI unit_2 = UnitSI().m().kg(0).s(-2);
	std::cout << unit_2.print() << std::endl;
	EXPECT_STREQ(std::string("$[ms^{-2}]$").c_str(), unit_2.print().c_str());

	UnitSI unit_3 = UnitSI().m(0); // dimensionless quantity
	std::cout << unit_3.print() << std::endl;
	EXPECT_STREQ(std::string("$[-]$").c_str(), unit_3.print().c_str());
}

TEST(UnitSI, multiplicative_operator) {
	UnitSI pressure = UnitSI().m(-1).kg().s(-2);
	UnitSI area = UnitSI().m(2);
	UnitSI force = pressure * area;
	EXPECT_STREQ(std::string("$[mkgs^{-2}]$").c_str(), force.print().c_str());

	UnitSI force2 = area * UnitSI::Pa;
	EXPECT_STREQ( UnitSI::N.print().c_str(), force2.print().c_str() );
}
