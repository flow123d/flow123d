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
#include "fields/unit_converter.hh"

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


template< class Iter_type >
class Semantic_actions
{
public:
    Semantic_actions()
    : unit_data_key_(""), factor_idx_(-1)
    {
    	unit_data_[""] = Formula();
    }

    void new_shortcut( Iter_type begin, Iter_type end )
    {
    	std::string key = get_str(begin, end);
    	unit_data_[key] = Formula();
    	unit_data_key_ = key;
    	factor_idx_ = -1;
    }

    void new_mult_factor( Iter_type begin, Iter_type end )
    {
    	unit_data_[unit_data_key_].factors_.push_back( Factor(get_str(begin, end), 1 ) );
    	++factor_idx_;
    }

    void new_div_factor( Iter_type begin, Iter_type end )
    {
    	unit_data_[unit_data_key_].factors_.push_back( Factor(get_str(begin, end), -1 ) );
    	++factor_idx_;
    }

    void new_exp( boost::int64_t i )
    {
    	unit_data_[unit_data_key_].factors_[factor_idx_].exponent_ *= (int)i;
    }

    void new_multipl( double d )
    {
    	unit_data_[unit_data_key_].coef_ = d;
    }
private:

    Semantic_actions& operator=( const Semantic_actions& );
                                // to prevent "assignment operator could not be generated" warning

    inline std::string get_str( std::string::const_iterator begin, std::string::const_iterator end )
    {
        return std::string( begin, end );
    }


    UnitData unit_data_;         //!< Full parsed data
    std::string unit_data_key_;  //!< keyo actual item of unit_data_
    int factor_idx_;             //!< index to actual item of subvector of factors of unit_data_
};


// the spirit grammer
template< class Iter_type >
class UnitSIGrammer : public spirit_namespace::grammar< UnitSIGrammer< Iter_type > >
{
public:

    typedef Semantic_actions< Iter_type > Semantic_actions_t;

    UnitSIGrammer( Semantic_actions_t& semantic_actions )
    :   actions_( semantic_actions )
    {
    }

	static void throw_not_equating( Iter_type begin, Iter_type end )
    {
	    throw_error( begin, "not expression, missing '='" );
    }

	static void throw_exp_not_int( Iter_type begin, Iter_type end )
    {
	    throw_error( begin, "not integer value of exponent" );
    }


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

            Str_action    new_shortcut     ( boost::bind( &Semantic_actions_t::new_shortcut,     &self.actions_, _1, _2 ) );
            Str_action    new_mult_factor  ( boost::bind( &Semantic_actions_t::new_mult_factor,  &self.actions_, _1, _2 ) );
            Str_action    new_div_factor   ( boost::bind( &Semantic_actions_t::new_div_factor,   &self.actions_, _1, _2 ) );
            Real_action   new_multipl      ( boost::bind( &Semantic_actions_t::new_multipl,      &self.actions_, _1 ) );
            Int_action    new_exp          ( boost::bind( &Semantic_actions_t::new_exp,          &self.actions_, _1 ) );

            // actual grammer

            unit_
			    = formula_ >> *( ( ch_p(';') >> *space_p ) >> constant_ )
                ;

            formula_
                = formula_factor_[ new_mult_factor ] >> !( ch_p('^') >> exp_ )
                  >> *( (ch_p('*') >> formula_factor_[ new_mult_factor ] >> !( ch_p('^') >> exp_ ) )
                      | (ch_p('/') >> formula_factor_[ new_div_factor ] >> !( ch_p('^') >> exp_ ) )
                      )
                ;

            formula_factor_
		        = +alpha_p
			    ;

            constant_
			    = constant_shortcut_[ new_shortcut ]
				  >> *space_p >> ( ch_p('=') | eps_p[ &throw_not_equating ] )
				  >> *space_p >> !(constant_multiplicator_ >> ch_p('*'))
			      >> formula_
			    ;

            constant_shortcut_
			    = +alpha_p
				;

            constant_multiplicator_
			    = strict_real_p[ new_multipl ]
			    ;

            exp_
			    = int64_p[ new_exp ]
			      | eps_p[ &throw_exp_not_int ]
				;
        }

        spirit_namespace::rule< ScannerT > unit_, formula_, formula_factor_, exp_, constant_, constant_shortcut_, constant_multiplicator_;

        const spirit_namespace::rule< ScannerT >& start() const { return unit_; }
    };

private:

    UnitSIGrammer& operator=( const UnitSIGrammer& ); // to prevent "assignment operator could not be generated" warning

    Semantic_actions_t& actions_;
};

} // end namespace units_converter

#endif /* UNIT_CONVERTER_TEMPLATE_HH_ */
