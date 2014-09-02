/*
 * unit_si_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>

#include "fields/unit_si.hh"


TEST(UnitSI, format) {
	UnitSI unit_1 = UnitSI().m(2).kg(1).s(-2);
	EXPECT_STREQ(std::string("$[m^{2}kgs^{-2}]$").c_str(), unit_1.format().c_str());

	UnitSI unit_2 = UnitSI().m().kg(0).s(-2);
	EXPECT_STREQ(std::string("$[ms^{-2}]$").c_str(), unit_2.format().c_str());

	UnitSI unit_3 = UnitSI().m(0); // dimensionless quantity
	EXPECT_STREQ(std::string("$[-]$").c_str(), unit_3.format().c_str());
}

TEST(UnitSI, static_defined_units) {
	EXPECT_STREQ( std::string("$[mkgs^{-2}]$").c_str(), UnitSI::N().format().c_str() );
	EXPECT_STREQ( std::string("$[m^{2}kgs^{-2}]$").c_str(), UnitSI::J().format().c_str() );
	EXPECT_STREQ( std::string("$[m^{2}kgs^{-3}]$").c_str(), UnitSI::W().format().c_str() );
	EXPECT_STREQ( std::string("$[m^{-1}kgs^{-2}]$").c_str(), UnitSI::Pa().format().c_str() );
	EXPECT_STREQ( std::string("$[-]$").c_str(), UnitSI::dimensionless().format().c_str() );
}

TEST(UnitSI, multiplicative_operator) {
	UnitSI pressure = UnitSI().m(-1).kg().s(-2);
	UnitSI area = UnitSI().m(2);
	UnitSI force = pressure * area;
	EXPECT_STREQ(std::string("$[mkgs^{-2}]$").c_str(), force.format().c_str());

	UnitSI force2 = area * UnitSI::Pa();
	EXPECT_STREQ( UnitSI::N().format().c_str(), force2.format().c_str() );
}

TEST(UnitSI, division_operator) {
	UnitSI force = UnitSI().m().kg().s(-2);
	UnitSI area = UnitSI().m(2);
	UnitSI pressure = force / area;
	EXPECT_STREQ(std::string("$[m^{-1}kgs^{-2}]$").c_str(), pressure.format().c_str());

	UnitSI pressure2 = UnitSI::N() / area;
	EXPECT_STREQ( UnitSI::Pa().format().c_str(), pressure2.format().c_str() );
}
