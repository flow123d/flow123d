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
 * @file    isotherm.cc
 * @brief   
 */

#include "reaction/isotherm.hh"
#include "system/sys_profiler.hh"
#include "system/logger.hh"

#include <boost/core/explicit_operator_bool.hpp>              // for optiona...
#include <boost/exception/detail/error_info_impl.hpp>         // for error_info
#include <boost/exception/exception.hpp>                      // for exception
#include <boost/exception/info.hpp>                           // for error_i...
#include <boost/format/alt_sstream.hpp>                       // for basic_a...
#include <boost/format/alt_sstream_impl.hpp>                  // for basic_a...
#include <boost/optional/optional.hpp>                        // for get_poi...
#include <boost/exception/diagnostic_information.hpp>         // for diagnos...
#include <boost/math/tools/toms748_solve.hpp>                 // for toms748...

Isotherm::Isotherm()
: table_limit_(0.0)
{
}

void Isotherm::clear_table()
{
    table_limit_ = 0.0;
    interpolation_table.clear();
}


void Isotherm::reinit(enum SorptionType adsorption_type, bool limited_solubility_on,
                      double rho_aqua, double scale_aqua, double scale_sorbed,
                      double c_aqua_limit, double mult_coef, double second_coef)
{
    adsorption_type_ = adsorption_type;
    rho_aqua_ = rho_aqua;
    scale_aqua_ = scale_aqua;
    scale_sorbed_ = scale_sorbed;
    limited_solubility_on_ = limited_solubility_on;
    mult_coef_ = mult_coef*rho_aqua;
    second_coef_ = second_coef;
    
    solubility_limit_ = c_aqua_limit;
    inv_scale_aqua_ = scale_aqua_/(scale_aqua_*scale_aqua_ + scale_sorbed_*scale_sorbed_);
    inv_scale_sorbed_ = scale_sorbed_/(scale_aqua_*scale_aqua_ + scale_sorbed_*scale_sorbed_);
}


double Isotherm::get_total_mass( ConcPair conc )
{
    return scale_aqua_* conc.fluid + scale_sorbed_ * conc.solid;
}


void Isotherm::compute( double &c_aqua, double &c_sorbed )
{
    // if sorption is switched off, do not compute anything
    if(adsorption_type_ == SorptionType::none)
        return;
    
    ConcPair c_pair(c_aqua, c_sorbed);
    ConcPair result(0,0);

    result = solve_conc( c_pair );

    c_aqua=result.fluid;
    c_sorbed=result.solid;
}


void Isotherm::interpolate( double &c_aqua, double &c_sorbed )
{
    // if sorption is switched off, do not compute anything
    if(adsorption_type_ == SorptionType::none)
        return;
    
    ConcPair c_pair(c_aqua, c_sorbed);
    ConcPair result(0,0);

    result = compute_projection( c_pair );

    c_aqua=result.fluid;
    c_sorbed=result.solid;
}


Isotherm::ConcPair Isotherm::compute_projection( Isotherm::ConcPair c_pair )
{
    double total_mass = get_total_mass(c_pair);
//     DebugOut().fmt("compute_projection: total mass = {}, c_aqua = {}\n", total_mass, c_pair.fluid);
    if(total_mass < 0.0)
        THROW( Isotherm::ExcNegativeTotalMass() 
                << EI_TotalMass(total_mass)
                );
    // total_mass_step_ is set and checked in make_table
    double total_mass_steps = total_mass / total_mass_step_;
    unsigned int total_mass_idx = static_cast <unsigned int>(std::floor(total_mass_steps));

    if (total_mass_idx < (interpolation_table.size() - 1) ) {
        double rot_sorbed = interpolation_table[total_mass_idx]
                            + (total_mass_steps - total_mass_idx)*(interpolation_table[total_mass_idx+1]
                            - interpolation_table[total_mass_idx]);
        return ConcPair( (total_mass * inv_scale_aqua_ - rot_sorbed * inv_scale_sorbed_),
                         (total_mass * inv_scale_sorbed_ + rot_sorbed * inv_scale_aqua_) );
    } else {
        if (limited_solubility_on_)
            return precipitate( c_pair );
        else
            return solve_conc( c_pair );
  }
}


Isotherm::ConcPair Isotherm::precipitate( Isotherm::ConcPair c_pair)
{
    double total_mass = get_total_mass(c_pair);
//     DebugOut().fmt("precipitate: total mass = {}, c_aqua = {}\n", total_mass, c_pair.fluid);
    return ConcPair(solubility_limit_,
                    (total_mass - scale_aqua_ * solubility_limit_) / scale_sorbed_  );
}


template<class Func>
Isotherm::ConcPair Isotherm::solve_conc( Isotherm::ConcPair c_pair, const Func &isotherm )
{
    double total_mass = get_total_mass(c_pair);
    double mass_limit = get_total_mass(Isotherm::ConcPair(solubility_limit_, const_cast<Func &>(isotherm)(solubility_limit_ / this->rho_aqua_)));
    
    // condition on limited solubility in the rotated coordinate system (total mass)
//     DebugOut().fmt("total_mass {}, mass_limit {} \n",total_mass, mass_limit);
    if (total_mass > mass_limit){
        if(limited_solubility_on_)
            return precipitate( c_pair );
        else
            // if solubility is not limited, increase mass limit
            mass_limit = total_mass;
    }

    double upper_solution_bound = mass_limit / scale_aqua_;
    CrossFunction<Func> eq_func(isotherm, total_mass, scale_aqua_, scale_sorbed_, this->rho_aqua_);
    std::pair<double,double> solution;
    if (total_mass > 0) // here should be probably some kind of tolerance instead of "0"
    {
        try {
            boost::uintmax_t max_iter = 20;
            tolerance<double> toler(30);
            solution = boost::math::tools::toms748_solve(eq_func, 0.0, upper_solution_bound, toler, max_iter);
        }
        catch(boost::exception const & e)
        {
            THROW( Isotherm::ExcBoostSolver() 
                << EI_BoostMessage(boost::diagnostic_information(e))
                );
        }
    }
    double difference;
    difference = (solution.second - solution.first)/2;
    double c_aqua = solution.first + difference;
    return ConcPair(c_aqua,
                    (total_mass - scale_aqua_ * c_aqua)/scale_sorbed_);
}

// Isotherm None specialization
template<> Isotherm::ConcPair Isotherm::solve_conc( Isotherm::ConcPair c_pair, const None &)
{
    return c_pair;
}


Isotherm::ConcPair Isotherm::solve_conc( Isotherm::ConcPair conc )
{
    switch(adsorption_type_)
    {
        case 0: // none
//             {
//                     None obj_isotherm;
//                     return solve_conc( conc, obj_isotherm);
//             }
            break;
        case 1: //  linear:
            {
                    Linear obj_isotherm(mult_coef_);
                    return solve_conc( conc, obj_isotherm);
            }
            break;
        case 2: // freundlich
            {
                    Freundlich obj_isotherm(mult_coef_, second_coef_);
                    return solve_conc( conc, obj_isotherm);
            }
            break;
        case 3:  // langmuir:
            {
                    Langmuir obj_isotherm(mult_coef_, second_coef_);
                    return solve_conc( conc, obj_isotherm);
            }
            break;
    }
    return conc;
}


template<class Func>
void Isotherm::make_table( const Func &isotherm, int n_steps )
{
    // limit aqueous concentration for the interpolation table; cannot be higher than solubility limit
    double aqua_limit = table_limit_;
    if(limited_solubility_on_)
            aqua_limit = solubility_limit_;
    
    double mass_limit = scale_aqua_ * aqua_limit + scale_sorbed_ * const_cast<Func &>(isotherm)(aqua_limit / this->rho_aqua_);
//     DebugOut().fmt("make_table: mass_limit = {}, aqua_limit = {}\n", mass_limit, aqua_limit);
    
    if(mass_limit < 0.0)
        THROW( Isotherm::ExcNegativeTotalMass()
                << EI_TotalMass(mass_limit)
                );
    total_mass_step_ = mass_limit / n_steps;
    double mass = 0.0;
    interpolation_table.clear();
    interpolation_table.reserve(n_steps);
    for(int i=0; i<= n_steps; i++) {
         // aqueous concentration (original coordinates c_a) corresponding to i-th total_mass_step_
        ConcPair c_pair( mass/scale_aqua_, 0.0 );

        ConcPair result = solve_conc( c_pair, isotherm);
        double c_sorbed_rot = ( result.solid * scale_aqua_ - result.fluid * scale_sorbed_);
        interpolation_table.push_back(c_sorbed_rot);
        mass = mass+total_mass_step_;
    }
}

// Isotherm None specialization
template<> void Isotherm::make_table( const None &, int )
{
    // Solve_conc returns the same, so we need to do that also in compute_projection.
    // We set size of the table to 1, so it follow the conditions into solve_conc again.
    
    limited_solubility_on_ = false; // so it cannot go in precipitate function
    
    total_mass_step_ = 1;            // set just one step in the table, so we void zero division
    interpolation_table.resize(1,0); // set one value in the table so the condition in compute_projection fails
    return;
}

void Isotherm::make_table( unsigned int n_points, double table_limit )
{
    START_TIMER("Isotherm::make_table");
    table_limit_ = table_limit;
    if(table_limit_ > 0.0) 
        switch(adsorption_type_)
    {
        case 0: // none
            {
                    None obj_isotherm;
                    make_table(obj_isotherm, 1);
            }
            break;
        case 1: //  linear:
            {
                Linear obj_isotherm(mult_coef_);
                make_table(obj_isotherm, n_points);
            }
            break;
        case 2: // freundlich:
            {
                Freundlich obj_isotherm(mult_coef_, second_coef_);
                make_table(obj_isotherm, n_points);
            }
            break;
        case 3: // langmuir:
            {
                Langmuir obj_isotherm(mult_coef_, second_coef_);
                make_table(obj_isotherm, n_points);
            }
            break;
        default:
            break;
    }
    else
        clear_table();
}
