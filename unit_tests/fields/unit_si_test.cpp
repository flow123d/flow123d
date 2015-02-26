/*
 * unit_si_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>

#include "fields/unit_si.hh"


TEST(UnitSI, format_latex) {
	UnitSI unit_1 = UnitSI().m(2).kg(1).s(-2);
	EXPECT_EQ("m^{2}kgs^{-2}", unit_1.format_latex());
    EXPECT_EQ("m(2).kg.s(-2)", unit_1.format_text());


	UnitSI unit_2 = UnitSI().m().kg(0).s(-2);
	EXPECT_EQ("ms^{-2}", unit_2.format_latex());
	EXPECT_EQ("m.s(-2)", unit_2.format_text());

	UnitSI unit_3 = UnitSI().m(0); // dimensionless quantity
	EXPECT_EQ("-", unit_3.format_latex());
	EXPECT_EQ("-", unit_3.format_text());

	UnitSI unit_4 = UnitSI().m(4).md();
	EXPECT_EQ("m^{4-d}", unit_4.format_latex());
	EXPECT_EQ("m(4-d)", unit_4.format_text());

	UnitSI unit_5 = UnitSI().m().md(1);
	EXPECT_EQ("m^{1+d}", unit_5.format_latex());
	EXPECT_EQ("m(1+d)", unit_5.format_text());

	UnitSI unit_7 = UnitSI().m(-1).md(1);
    EXPECT_EQ("m^{-1+d}", unit_7.format_latex());
    EXPECT_EQ("m(-1+d)", unit_7.format_text());

	UnitSI unit_6 = UnitSI().kg().md(1);
	EXPECT_EQ("m^{d}kg", unit_6.format_latex());
	EXPECT_EQ("m(d).kg", unit_6.format_text());
}

TEST(UnitSI, static_defined_units) {
	EXPECT_EQ("mkgs^{-2}", UnitSI::N().format_latex());
	EXPECT_EQ("m^{2}kgs^{-2}", UnitSI::J().format_latex());
	EXPECT_EQ("m^{2}kgs^{-3}", UnitSI::W().format_latex());
	EXPECT_EQ("m^{-1}kgs^{-2}", UnitSI::Pa().format_latex());
	EXPECT_EQ("-", UnitSI::dimensionless().format_latex());
}

TEST(UnitSI, multiplicative_operator) {
	UnitSI pressure = UnitSI().m(-1).kg().s(-2);
	UnitSI area = UnitSI().m(2);
	UnitSI force = pressure * area;
	EXPECT_EQ("mkgs^{-2}", force.format_latex());

	UnitSI force2 = area * UnitSI::Pa();
	EXPECT_EQ( UnitSI::N().format_latex(), force2.format_latex() );
}

TEST(UnitSI, division_operator) {
	UnitSI force = UnitSI().m().kg().s(-2);
	UnitSI area = UnitSI().m(2);
	UnitSI pressure = force / area;
	EXPECT_EQ("m^{-1}kgs^{-2}", pressure.format_latex());

	UnitSI pressure2 = UnitSI::N() / area;
	EXPECT_EQ( UnitSI::Pa().format_latex(), pressure2.format_latex() );
}
