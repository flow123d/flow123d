/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 *
 * @file    unit_converter_template.hh
 * @brief
 */

#ifndef UNIT_CONVERTER_TEMPLATE_HH_
#define UNIT_CONVERTER_TEMPLATE_HH_


#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/version.hpp>
#include "tools/unit_converter.hh"

#if BOOST_VERSION >= 103800
    #include <boost/spirit/include/classic_core.hpp>
    #include <boost/spirit/include/classic_confix.hpp>
    #include <boost/spirit/include/classic_escape_char.hpp>
    #include <boost/spirit/include/classic_multi_pass.hpp>
    #include <boost/spirit/include/classic_position_iterator.hpp>
    #define spirit_namespace boost::spirit::classic
#else
    #include <boost/spirit/core.hpp>
    #include <boost/spirit/utility/confix.hpp>
    #include <boost/spirit/utility/escape_char.hpp>
    #include <boost/spirit/iterator/multi_pass.hpp>
    #include <boost/spirit/iterator/position_iterator.hpp>
    #define spirit_namespace boost::spirit
#endif


namespace units_converter
{

const spirit_namespace::int_parser < boost::int64_t >  int64_p  = spirit_namespace::int_parser < boost::int64_t  >();
const spirit_namespace::uint_parser< boost::uint64_t > uint64_p = spirit_namespace::uint_parser< boost::uint64_t >();


/**
 * @brief Class manages parsing of user defined field unit.
 *
 * Class contains:
 *  - object of type \p UnitData for storing user defined unit
 *  - methods for create UnitData object
 *  - methods managing exceptions, if user defined unit is not in correct format
 *
 * For example, unit is defined in format:
 * "MPa/rho/g_; rho = 990*kg*m^-3; g_ = 9.8*m*s^-2"
 *
 *  - this unit is composed of three parts separated by semicolons
 *  - first part (MPa/rho/g) defines unit and allows operations multiplication, division and exponentiation of factors
 *  - factor 'MPa' is predefined (MegaPascal), factors 'rho' and 'g_' must be defined in next parts
 *  - first part is always obligatory
 *  - following parts must be in format: '<shortcut> = (<multipicative_coeficient>*)<definition>'
 *  - shotcut must correspond with user defined factor in first part
 *  - multipicative_coeficient is optional
 *  - definition allows operations multiplication, division and exponentiation of factors too
 *  - symbol 'g_' must be defined in this format because 'g' is in conflict with gram
 */
template< class Iter_type >
class Semantic_actions
{
public:
	/// Constructor.
    Semantic_actions()
    : unit_data_key_(""), factor_idx_(-1)
    {
    	unit_data_[""] = Formula();
    }

    /// Add new definition of formula
    void new_shortcut( Iter_type begin, Iter_type end )
    {
    	std::string key = get_str(begin, end);
    	unit_data_[key] = Formula();
    	unit_data_key_ = key;
    	factor_idx_ = -1;
    }

    /// Add new factor of unit (factor is multiplying)
    void new_mult_factor( Iter_type begin, Iter_type end )
    {
    	unit_data_[unit_data_key_].factors_.push_back( Factor(get_str(begin, end), 1 ) );
    	++factor_idx_;
    }

    /// Add new factor of unit (factor is dividing)
    void new_div_factor( Iter_type begin, Iter_type end )
    {
    	unit_data_[unit_data_key_].factors_.push_back( Factor(get_str(begin, end), -1 ) );
    	++factor_idx_;
    }

    /// Compute exponent to actual factor of unit
    void new_exp( boost::int64_t i )
    {
    	unit_data_[unit_data_key_].factors_[factor_idx_].exponent_ *= (int)i;
    }

    /// Add multipicative coeficient of unit
    void new_multipl( double d )
    {
    	unit_data_[unit_data_key_].coef_ = d;
    }

    /// Throw exception if exponent is not in correct format
	void throw_exp_not_int( Iter_type begin, Iter_type end )
    {
		std::stringstream ss;
		ss << "Value of exponent '" << get_str( begin, end ) << "' is not integer";
		THROW( ExcInvalidUnit() << EI_UnitError(ss.str()) );
    }

	/// Throw exception if shortcut of factor is not in correct format
	void throw_not_shortcut( Iter_type begin, Iter_type end )
    {
		std::stringstream ss;
		if (begin == end) {
			ss << "Missing declaration of shortcut '.." << get_str( begin-4, end+4 ) << "..'";
		} else {
			ss << "Invalid shortcut of unit '" << get_str( begin, end ) << "'";
		}
		THROW( ExcInvalidUnit() << EI_UnitError(ss.str()) );
    }

	/// Throw exception if sign '=' missing in definition
	void throw_not_equating( Iter_type begin, Iter_type end )
    {
		std::stringstream ss;
		ss << "Invalid expression '" << get_str( begin, end ) << "', missing '='";
		THROW( ExcInvalidUnit() << EI_UnitError(ss.str()) );
    }

	/**
	 * @brief Check @p unit_data_ object.
	 *
	 * Method:
	 *  - marks factors that are defined as derived unit
	 *  - checks undefined factors of unit
	 *  - check conflicts in definitions of unit (same shortcut is defined by user and predefined in application)
	 *  - checks cyclic definition of unit
	 */
    void check_unit_data()
    {
    	for(std::map<std::string, struct Formula>::iterator it = unit_data_.begin(); it != unit_data_.end(); ++it)
    	{
    		for(std::vector<struct Factor>::iterator formula_it = it->second.factors_.begin();
    				formula_it != it->second.factors_.end(); ++formula_it)
    		{
    			if (formula_it->factor_ == it->first) { // error - cyclic definition
    				std::stringstream ss;
    				ss << "Cyclic declaration of unit '" << it->first << "'";
    				THROW( ExcInvalidUnit() << EI_UnitError(ss.str()) );
    			} else if (unit_data_.find(formula_it->factor_) != unit_data_.end()) { // formula exists as derived unit
    				if (UnitConverter::basic_factors.units_map_.find(formula_it->factor_) != UnitConverter::basic_factors.units_map_.end()) {
    					// error - conflict with predefined unit
        				std::stringstream ss;
        				ss << "Shortcut '" << formula_it->factor_ << "' is in conflict with predefined unit";
        				THROW( ExcInvalidUnit() << EI_UnitError(ss.str()) );
    				}
    				formula_it->basic_ = false;
    			} else if (UnitConverter::basic_factors.units_map_.find(formula_it->factor_) == UnitConverter::basic_factors.units_map_.end()) {
    				// error - unit not defined
    				std::stringstream ss;
    				ss << "Unit '" << formula_it->factor_ << "' is not defined";
    				THROW( ExcInvalidUnit() << EI_UnitError(ss.str()) );
    			}
    		}
    	}
    }

    // Return @p unit_data_
    inline UnitData unit_data() const
    { return unit_data_; }
private:

    Semantic_actions& operator=( const Semantic_actions& );
                                // to prevent "assignment operator could not be generated" warning

    inline std::string get_str( std::string::const_iterator begin, std::string::const_iterator end ) const
    {
        return std::string( begin, end );
    }


    UnitData unit_data_;         //!< Full parsed data
    std::string unit_data_key_;  //!< key of actual item of unit_data_
    int factor_idx_;             //!< index to actual item of subvector of factors of unit_data_
};


/**
 * @brief Definition of unit grammar.
 *
 * Allow parse user-defined units.
 */
template< class Iter_type >
class UnitSIGrammer : public spirit_namespace::grammar< UnitSIGrammer< Iter_type > >
{
public:

    typedef Semantic_actions< Iter_type > Semantic_actions_t;

    /// Constructor.
    UnitSIGrammer( Semantic_actions_t& semantic_actions )
    :   actions_( semantic_actions )
    {
    }


    /// Define rules of grammar
    template< typename ScannerT >
    class definition
    {
    public:

        definition( const UnitSIGrammer& self )
        {
            using namespace spirit_namespace;

            // first we convert the semantic action class methods to functors with the
            // parameter signature expected by spirit

            typedef boost::function< void( Iter_type, Iter_type ) > Str_action;
            typedef boost::function< void( double )               > Real_action;
            typedef boost::function< void( boost::int64_t )       > Int_action;

            Str_action    new_shortcut       ( boost::bind( &Semantic_actions_t::new_shortcut,       &self.actions_, _1, _2 ) );
            Str_action    new_mult_factor    ( boost::bind( &Semantic_actions_t::new_mult_factor,    &self.actions_, _1, _2 ) );
            Str_action    new_div_factor     ( boost::bind( &Semantic_actions_t::new_div_factor,     &self.actions_, _1, _2 ) );
            Real_action   new_multipl        ( boost::bind( &Semantic_actions_t::new_multipl,        &self.actions_, _1 ) );
            Int_action    new_exp            ( boost::bind( &Semantic_actions_t::new_exp,            &self.actions_, _1 ) );
            Str_action    throw_exp_not_int  ( boost::bind( &Semantic_actions_t::throw_exp_not_int,  &self.actions_, _1, _2 ) );
            Str_action    throw_not_equating ( boost::bind( &Semantic_actions_t::throw_not_equating, &self.actions_, _1, _2 ) );
            Str_action    throw_not_shortcut ( boost::bind( &Semantic_actions_t::throw_not_shortcut, &self.actions_, _1, _2 ) );

            // actual grammer

            unit_
			    = formula_ >> *( ch_p(';') >> *space_p >> constant_ )
                ;

            formula_
                = formula_factor_[ new_mult_factor ] >> !( ch_p('^') >> exp_ )
                  >> *( (ch_p('*') >> formula_factor_[ new_mult_factor ] >> !( ch_p('^') >> exp_ ) )
                      | (ch_p('/') >> formula_factor_[ new_div_factor ] >> !( ch_p('^') >> exp_ ) )
                      )
                ;

            formula_factor_
		        = forbidden_char_[ throw_not_shortcut ]
		          | ( alpha_p >> *( anychar_p - ch_p(';') - ch_p('*') - ch_p('/') - ch_p('=') - ch_p('^') - space_p ) )
		          | shortcut_err_[ throw_not_shortcut ]
			    ;

            constant_
			    = ( formula_factor_[ new_shortcut ] >> *space_p
			        >> ch_p('=') >> *space_p >> constant_value_ )
				  | constant_err_[ throw_not_equating ]
			    ;

            constant_value_
			    = !(constant_multiplicator_ >> ch_p('*'))
				  >> formula_
				;

            constant_err_
		        = +( anychar_p - ch_p(';') )
				;

            constant_multiplicator_
			    = strict_real_p[ new_multipl ]
				  | int64_p[ new_multipl ]
			    ;

            exp_
			    = int64_p[ new_exp ]
			      | exp_err_[ throw_exp_not_int ]
				;

            exp_err_
			    = +( anychar_p - ch_p(';') - ch_p('*') - ch_p('/') - space_p )
				;

            shortcut_err_
			    = *( anychar_p - ch_p(';') - ch_p('*') - ch_p('/') - ch_p(';') - ch_p('^') )
				;

            forbidden_char_
			    = *( alnum_p | ch_p('_') )
				>> ( anychar_p - alnum_p - ch_p('_') - ch_p(';') - ch_p('*') - ch_p('/') - ch_p('=') - ch_p('^') - space_p )
				>> *( alnum_p | ch_p('_') )
                ;

        }

        spirit_namespace::rule< ScannerT > unit_, formula_, formula_factor_, exp_, exp_err_, constant_,
				constant_value_, constant_err_, constant_multiplicator_, shortcut_err_, forbidden_char_;

        const spirit_namespace::rule< ScannerT >& start() const { return unit_; }
    };

private:

    UnitSIGrammer& operator=( const UnitSIGrammer& ); // to prevent "assignment operator could not be generated" warning

    Semantic_actions_t& actions_;
};


} // end namespace units_converter

#endif /* UNIT_CONVERTER_TEMPLATE_HH_ */
