/*
 * unit_si_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include "tools/unit_si.hh"
#include "tools/unit_converter.hh"
#include "tools/unit_converter_template.hh"


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


class UnitConverterTest : public testing::Test, public UnitConverter {
public:
	UnitData unit_data_;
protected:

    virtual void SetUp() {
    }
    virtual void TearDown() {
    };

};


TEST_F(UnitConverterTest, converter_grammar) {
	{
		std::string unit = "MPa/rho/g_; rho = 990*kg*m^-3; g_ = 9.8*m*s^-2";
		unit_data_ = this->read_unit(unit);
		EXPECT_EQ(unit_data_.size(), 3);
		EXPECT_TRUE( unit_data_.find("rho") != unit_data_.end() );
		EXPECT_TRUE( unit_data_.find("g_") != unit_data_.end() );
		EXPECT_FALSE( unit_data_.find("MPa") != unit_data_.end() );
		EXPECT_DOUBLE_EQ( unit_data_.find("g_")->second.coef_, 9.8);
		EXPECT_EQ(unit_data_.find("g_")->second.factors_.size(), 2);
		EXPECT_DOUBLE_EQ( unit_data_.find("rho")->second.coef_, 990);
		EXPECT_EQ(unit_data_.find("rho")->second.factors_.size(), 2);
	}

	{
		std::string unit = "a*b; a = kg*m^a; b = m*s^-2";
		EXPECT_THROW_WHAT( { unit_data_ = this->read_unit(unit); }, ExcInvalidUnit,
				"Value of exponent 'a' is not integer" );
	}

	{
		std::string unit = "a*b; b*c; b = m*s^-2";
		EXPECT_THROW_WHAT( { unit_data_ = this->read_unit(unit); }, ExcInvalidUnit,
				"Invalid expression '.*', missing '='" );
	}

	{
		std::string unit = "kPa*n; n = m*55";
		EXPECT_THROW_WHAT( { unit_data_ = this->read_unit(unit); }, ExcInvalidUnit,
				"Invalid shortcut of unit '55'" );
	}

	{
		std::string unit = "a*b; a = m*; b = kg*s*m^2";
		EXPECT_THROW_WHAT( { unit_data_ = this->read_unit(unit); }, ExcInvalidUnit,
				"Missing declaration of shortcut '.*'" );
	}

	{
		std::string unit = "a*b; a = m*a; b = kg*s*m^2";
		EXPECT_THROW_WHAT( { unit_data_ = this->read_unit(unit); }, ExcInvalidUnit,
				"Cyclic declaration of unit 'a'" );
	}

	{
		std::string unit = "s*g; g = 9.8*m*s^-2";
		EXPECT_THROW_WHAT( { unit_data_ = this->read_unit(unit); }, ExcInvalidUnit,
				"Shortcut 'g' is in conflict with predefined unit" );
	}

	{
		std::string unit = "a*b*c; a = kg*K; b = m*s^-2";
		EXPECT_THROW_WHAT( { unit_data_ = this->read_unit(unit); }, ExcInvalidUnit,
				"Unit 'c' is not defined" );
	}

	{
		std::string unit = "a+b; a=m*s^-2; b=kg/m^2";
		EXPECT_THROW_WHAT( { unit_data_ = this->read_unit(unit); }, ExcInvalidUnit,
				"Invalid shortcut of unit '.*'" );
	}
}

TEST_F(UnitConverterTest, convert_function) {
	double coef;
	{
		std::string unit = "N/m^2";
		coef = this->convert(unit);
		EXPECT_DOUBLE_EQ( coef, 1.0 );
		EXPECT_TRUE( this->unit_si()==UnitSI::Pa() );
	}
	{
		std::string unit = "kN/cm^2";
		coef = this->convert(unit);
		EXPECT_DOUBLE_EQ( coef, 1.0e7 );
		EXPECT_TRUE( this->unit_si()==UnitSI::Pa() );
	}
	{
		std::string unit = "MPa/rho/g_; rho = 990*kg*m^-3; g_ = 9.8*m*s^-2";
		coef = this->convert(unit);
		EXPECT_DOUBLE_EQ( coef, (1000.0/0.99/9.8) );
		EXPECT_TRUE( this->unit_si()==UnitSI().m() );
	}
	{
		std::string unit = "g_^2; g_ = 9.8*m*s^-2";
		coef = this->convert(unit);
		EXPECT_DOUBLE_EQ( coef, (9.8*9.8) );
		EXPECT_TRUE( this->unit_si()==UnitSI().m(2).s(-4) );
	}
	{
		std::string unit = "V*rho; rho = g*cm^-3; V = m^3";
		coef = this->convert(unit);
		EXPECT_DOUBLE_EQ( coef, 1000.0 );
		EXPECT_TRUE( this->unit_si()==UnitSI().kg() );
	}
	{
		std::string unit = "q/rho^2; q = 100*kg ; rho = 1000*kg*m^-3";
		coef = this->convert(unit);
		EXPECT_DOUBLE_EQ( coef, 0.0001 );
		EXPECT_TRUE( this->unit_si()==UnitSI().kg(-1).m(6) );
	}
}


TEST(UnitSI, convert_from_string) {
	{
		UnitSI unit = UnitSI().m();
		double coef = unit.convert_unit_from("MPa/rho/g_; rho = 990*kg*m^-3; g_ = 9.8*m*s^-2");
		EXPECT_DOUBLE_EQ( coef, (1000.0/0.99/9.8) );
	}
	{
		UnitSI unit = UnitSI().kg().m(-3);
		double coef = unit.convert_unit_from("g*cm^-3");
		EXPECT_DOUBLE_EQ( coef, 1000.0 );
	}
	{
		UnitSI unit = UnitSI().kg().m(-3);
		EXPECT_THROW_WHAT( { unit.convert_unit_from("m*rho; rho = g*cm^-3"); }, ExcNoncorrespondingUnit,
				"Non-corresponding definition of unit" );
	}
}
