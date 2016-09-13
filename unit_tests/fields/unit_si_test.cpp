/*
 * unit_si_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include "fields/unit_si.hh"
#include "fields/unit_converter_template.hh"


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

	UnitSI unit_8 = UnitSI().kg();
    EXPECT_EQ("kg", unit_8.format_latex());
    EXPECT_EQ("kg", unit_8.format_text());

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

TEST(UnitSI, user_defined_units) {
	{
		UnitSI unit = UnitSI("m");
		EXPECT_TRUE( unit == UnitSI().m() );
		EXPECT_DOUBLE_EQ(unit.coef(), 1);
	}

	{
		UnitSI unit = UnitSI("h");
		EXPECT_TRUE( unit == UnitSI().s() );
		EXPECT_DOUBLE_EQ(unit.coef(), 3600);
	}

	{
		UnitSI unit = UnitSI("kg.m^-3");
		EXPECT_TRUE( unit == UnitSI().m(-3).kg() );
		EXPECT_DOUBLE_EQ(unit.coef(), 1);
	}

	{
		UnitSI unit = UnitSI("g.cm^-3");
		EXPECT_TRUE( unit == UnitSI().m(-3).kg() );
		EXPECT_DOUBLE_EQ(unit.coef(), 1000);
	}

	{
		UnitSI unit = UnitSI("m.kg.s^-2");
		EXPECT_TRUE( unit == UnitSI().m().kg().s(-2) );
		EXPECT_DOUBLE_EQ(unit.coef(), 1);
	}

	{
		UnitSI unit = UnitSI("N.m^-2");
		EXPECT_TRUE( unit == UnitSI::Pa() );
		EXPECT_DOUBLE_EQ(unit.coef(), 1);
	}

	// Invalid representations of units
	EXPECT_THROW_WHAT( { UnitSI unit = UnitSI("m.s^-1^2"); }, UnitSI::ExcInvalidUnitString, "invalid value of unit" );
	EXPECT_THROW_WHAT( { UnitSI unit = UnitSI("ab^2"); }, UnitSI::ExcInvalidUnitString, "invalid symbol of unit 'ab'" );
	EXPECT_THROW_WHAT( { UnitSI unit = UnitSI("m^a"); }, UnitSI::ExcInvalidUnitString, "invalid exponent 'a'" );
	EXPECT_THROW_WHAT( { UnitSI unit = UnitSI("m^2a"); }, UnitSI::ExcInvalidUnitString, "invalid exponent '2a'" );
	EXPECT_THROW_WHAT( { UnitSI unit = UnitSI("m^1.5"); }, UnitSI::ExcInvalidUnitString, "invalid symbol of unit '5'" );
}

TEST(UnitSI, unit_converter_grammar) {
	units_converter::UnitData data;
	{
		std::string unit = "MPa/rho/g_; rho = 990*kg*m^-3; g_ = 9.8*m*s^-2";
		data = units_converter::read_unit(unit);
		EXPECT_EQ(data.size(), 3);
		EXPECT_TRUE( data.find("rho") != data.end() );
		EXPECT_TRUE( data.find("g_") != data.end() );
		EXPECT_FALSE( data.find("MPa") != data.end() );
		EXPECT_DOUBLE_EQ( data.find("g_")->second.coef_, 9.8);
		EXPECT_EQ(data.find("g_")->second.factors_.size(), 2);
		EXPECT_DOUBLE_EQ( data.find("rho")->second.coef_, 990);
		EXPECT_EQ(data.find("rho")->second.factors_.size(), 2);
	}

	{
		std::string unit = "a*b; a = kg*m^a; b = m*s^-2";
		EXPECT_THROW_WHAT( { data = units_converter::read_unit(unit); }, units_converter::ExcInvalidUnit,
				"Value of exponent 'a' is not integer" );
	}

	{
		std::string unit = "a*b; b*c; b = m*s^-2";
		EXPECT_THROW_WHAT( { data = units_converter::read_unit(unit); }, units_converter::ExcInvalidUnit,
				"Invalid expression '.*', missing '='" );
	}

	{
		std::string unit = "kPa*n; n = m*55";
		EXPECT_THROW_WHAT( { data = units_converter::read_unit(unit); }, units_converter::ExcInvalidUnit,
				"Invalid shortcut of unit '55'" );
	}

	{
		std::string unit = "a*b; a = m*; b = kg*s*m^2";
		EXPECT_THROW_WHAT( { data = units_converter::read_unit(unit); }, units_converter::ExcInvalidUnit,
				"Missing declaration of shortcut '.*'" );
	}

	{
		std::string unit = "a*b; a = m*a; b = kg*s*m^2";
		EXPECT_THROW_WHAT( { data = units_converter::read_unit(unit); }, units_converter::ExcInvalidUnit,
				"Cyclic declaration of unit 'a'" );
	}

	{
		std::string unit = "s*g; g = 9.8*m*s^-2";
		EXPECT_THROW_WHAT( { data = units_converter::read_unit(unit); }, units_converter::ExcInvalidUnit,
				"Shortcut 'g' is in conflict with predefined unit" );
	}

	{
		std::string unit = "a*b*c; a = kg.K; b = m*s^-2";
		EXPECT_THROW_WHAT( { data = units_converter::read_unit(unit); }, units_converter::ExcInvalidUnit,
				"Unit 'c' is not defined" );
	}
}
