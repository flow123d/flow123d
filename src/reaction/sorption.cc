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
 * @file    sorption.cc
 * @brief   
 */

#include <vector>
#include <limits>

#include "reaction/isotherm.hh"
#include "reaction/sorption.hh"
#include "system/sys_profiler.hh"
#include "system/asserts.hh"
#include "mesh/accessors.hh"
#include "input/factory.hh"

FLOW123D_FORCE_LINK_IN_CHILD(sorptionMobile)
FLOW123D_FORCE_LINK_IN_CHILD(sorptionImmobile)
FLOW123D_FORCE_LINK_IN_CHILD(sorption)


/********************************* SORPTION_SIMPLE *********************************************************/
/*********************************                 *********************************************************/

const IT::Record & SorptionSimple::get_input_type() {
	return IT::Record("Sorption", "Sorption model in the reaction term of transport.")
        .derive_from( ReactionTerm::it_abstract_term() )
        .copy_keys(SorptionBase::get_input_type())
        .declare_key("output", make_output_type("Sorption", "conc_solid", "Concentration solution in the solid phase."),
             IT::Default("{ \"fields\": [ \"conc_solid\" ] }"),
             "Setting of the fields output.")

		.close();
}

SorptionSimple::SorptionSimple(Mesh &init_mesh, Input::Record in_rec)
  : SorptionBase(init_mesh, in_rec)
{
	data_ = new EqData("conc_solid", "Concentration solution in the solid phase.");
    this->eq_data_ = data_;
}

const int SorptionSimple::registrar =
		Input::register_class< SorptionSimple, Mesh &, Input::Record >("Sorption") +
		SorptionSimple::get_input_type().size();

SorptionSimple::~SorptionSimple(void)
{}

void SorptionSimple::compute_common_ele_data(const ElementAccessor<3> &elem)
{
    double rock_density = data_->rock_density.value(elem.centre(),elem);
    double por_m = data_->porosity.value(elem.centre(),elem);
    
    this->common_ele_data.scale_aqua = por_m;
    this->common_ele_data.scale_sorbed = (1 - por_m) * rock_density;
    this->common_ele_data.no_sorbing_surface_cond = 1-por_m;
}


/***********************************               *********************************************************/
/*********************************** SORPTION_DUAL *********************************************************/
/***********************************               *********************************************************/

SorptionDual::SorptionDual(Mesh &init_mesh, Input::Record in_rec,
                           const string &output_conc_name,
                           const string &output_conc_desc)
    : SorptionBase(init_mesh, in_rec)
{
    data_ = new EqData(output_conc_name, output_conc_desc);
    *data_+=immob_porosity_
        .flags_add(FieldFlag::input_copy)
        .name("porosity_immobile")
		.set_limits(0.0);
    this->eq_data_ = data_;
}

SorptionDual::~SorptionDual(void)
{}

/**********************************                  *******************************************************/
/*********************************** SORPTION_MOBILE *******************************************************/
/**********************************                  *******************************************************/

const IT::Record & SorptionMob::get_input_type() {
	return IT::Record("SorptionMobile", "Sorption model in the mobile zone, following the dual porosity model.")
        .derive_from( ReactionTerm::it_abstract_mobile_term() )
        .copy_keys(SorptionBase::get_input_type())
        .declare_key("output", make_output_type("SorptionMobile", "conc_solid", "Concentration solution in the solid mobile phase."),
             IT::Default("{ \"fields\": [ \"conc_solid\" ] }"),
             "Setting of the fields output.")

		.close();
}


const int SorptionMob::registrar =
		Input::register_class< SorptionMob, Mesh &, Input::Record >("SorptionMobile") +
		SorptionMob::get_input_type().size();


SorptionMob::SorptionMob(Mesh &init_mesh, Input::Record in_rec)
    : SorptionDual(init_mesh, in_rec, "conc_solid", "Concentration solution in the solid mobile phase.")
{}


SorptionMob::~SorptionMob(void)
{}

void SorptionMob::compute_common_ele_data(const ElementAccessor<3> &elem)
{
    double rock_density = data_->rock_density.value(elem.centre(),elem);
    double por_m = data_->porosity.value(elem.centre(),elem);
    double por_imm = immob_porosity_.value(elem.centre(),elem);
    double phi = por_m/(por_m + por_imm);
    
    this->common_ele_data.scale_aqua = por_m;
    this->common_ele_data.scale_sorbed = phi * (1 - por_m - por_imm) * rock_density;
    this->common_ele_data.no_sorbing_surface_cond = 1-por_m-por_imm;
}


/***********************************                   *****************************************************/
/*********************************** SORPTION_IMMOBILE *****************************************************/
/***********************************                   *****************************************************/

const IT::Record & SorptionImmob::get_input_type() {
	return IT::Record("SorptionImmobile", "Sorption model in the immobile zone, following the dual porosity model.")
        .derive_from( ReactionTerm::it_abstract_immobile_term() )
        .copy_keys(SorptionBase::get_input_type())
        .declare_key("output", make_output_type("SorptionImmobile", "conc_immobile_solid", "Concentration solution in the solid immobile phase."),
             IT::Default("{ \"fields\": [ \"conc_immobile_solid\" ] }"),
             "Setting of the fields output.")

		.close();
}

const int SorptionImmob::registrar =
		Input::register_class< SorptionImmob, Mesh &, Input::Record >("SorptionImmobile") +
		SorptionImmob::get_input_type().size();

SorptionImmob::SorptionImmob(Mesh &init_mesh, Input::Record in_rec)
: SorptionDual(init_mesh, in_rec, "conc_immobile_solid", "Concentration solution in the solid immobile phase.")
{}

SorptionImmob::~SorptionImmob(void)
{}

void SorptionImmob::compute_common_ele_data(const ElementAccessor<3> &elem)
{
    double rock_density = data_->rock_density.value(elem.centre(),elem);
    double por_m = data_->porosity.value(elem.centre(),elem);
    double por_imm = immob_porosity_.value(elem.centre(),elem);
    double phi = por_m/(por_m + por_imm);
    
    this->common_ele_data.scale_aqua = por_m;
    this->common_ele_data.scale_sorbed = (1 - phi) * (1 - por_m - por_imm) * rock_density;
    this->common_ele_data.no_sorbing_surface_cond = 1-por_m-por_imm;
}
