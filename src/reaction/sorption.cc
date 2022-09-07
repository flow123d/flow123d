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
#include "fields/field_model.hh"

FLOW123D_FORCE_LINK_IN_CHILD(sorptionMobile)
FLOW123D_FORCE_LINK_IN_CHILD(sorptionImmobile)
FLOW123D_FORCE_LINK_IN_CHILD(sorption)


/********************************* SORPTION_SIMPLE *********************************************************/
/*********************************                 *********************************************************/

// Functors for computing models
struct fn_simple_scale_aqua {
    inline double operator() (double por_m) {
        return por_m;
    }
};

struct fn_simple_scale_sorbed {
    inline double operator() (double surface_cond, double rock_density) {
        return surface_cond * rock_density;
    }
};

struct fn_simple_surface_cond {
    inline double operator() (double por_m) {
        return 1 - por_m;
    }
};


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
    eq_fields_ = std::make_shared<EqFields>("conc_solid", "Concentration solution in the solid phase.");
    this->eq_fieldset_ = eq_fields_;
    this->eq_fields_base_ = std::static_pointer_cast<ReactionTerm::EqFields>(eq_fields_);
}

const int SorptionSimple::registrar =
		Input::register_class< SorptionSimple, Mesh &, Input::Record >("Sorption") +
		SorptionSimple::get_input_type().size();

SorptionSimple::~SorptionSimple(void)
{}

void SorptionSimple::init_field_models()
{
    eq_fields_->scale_aqua.set(Model<3, FieldValue<3>::Scalar>::create(fn_simple_scale_aqua(), eq_fields_->porosity), 0.0);
    eq_fields_->scale_sorbed.set(Model<3, FieldValue<3>::Scalar>::create(
            fn_simple_scale_sorbed(), eq_fields_->no_sorbing_surface_cond, eq_fields_->rock_density), 0.0);
    eq_fields_->no_sorbing_surface_cond.set(Model<3, FieldValue<3>::Scalar>::create(fn_simple_surface_cond(), eq_fields_->porosity), 0.0);
}


/***********************************               *********************************************************/
/*********************************** SORPTION_DUAL *********************************************************/
/***********************************               *********************************************************/

SorptionDual::EqFields::EqFields(const string &output_field_name, const string &output_field_desc)
: SorptionBase::EqFields(output_field_name, output_field_desc)
{
    *this+=immob_porosity_
            .flags_add(FieldFlag::input_copy)
            .name("porosity_immobile")
            .set_limits(0.0);
}

SorptionDual::SorptionDual(Mesh &init_mesh, Input::Record in_rec,
                           const string &output_conc_name,
                           const string &output_conc_desc)
    : SorptionBase(init_mesh, in_rec)
{
    eq_fields_dual_ = std::make_shared<EqFields>(output_conc_name, output_conc_desc);
    this->eq_fieldset_ = eq_fields_dual_;
    this->eq_fields_base_ = std::static_pointer_cast<ReactionTerm::EqFields>(eq_fields_dual_);
    this->eq_fields_ = std::static_pointer_cast<SorptionBase::EqFields>(eq_fields_dual_);
}

SorptionDual::~SorptionDual(void)
{}

/**********************************                  *******************************************************/
/*********************************** SORPTION_MOBILE *******************************************************/
/**********************************                  *******************************************************/

// Functors for computing models
struct fn_mob_scale_aqua {
    inline double operator() (double por_m) {
        return por_m;
    }
};

struct fn_mob_scale_sorbed {
    inline double operator() (double por_m, double por_imm, double surface_cond, double rock_density) {
        double phi = por_m/(por_m + por_imm);
    	return phi * surface_cond * rock_density;
    }
};

struct fn_mob_surface_cond {
    inline double operator() (double por_m, double por_imm) {
        return 1 - por_m - por_imm;
    }
};


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

void SorptionMob::init_field_models()
{
    eq_fields_->scale_aqua.set(Model<3, FieldValue<3>::Scalar>::create(fn_mob_scale_aqua(), eq_fields_->porosity), 0.0);
    eq_fields_->scale_sorbed.set(Model<3, FieldValue<3>::Scalar>::create(
            fn_mob_scale_sorbed(), eq_fields_->porosity, eq_fields_dual_->immob_porosity_, eq_fields_->no_sorbing_surface_cond,
	        eq_fields_->rock_density), 0.0);
    eq_fields_->no_sorbing_surface_cond.set(Model<3, FieldValue<3>::Scalar>::create(
            fn_mob_surface_cond(), eq_fields_->porosity, eq_fields_dual_->immob_porosity_), 0.0);
}


/***********************************                   *****************************************************/
/*********************************** SORPTION_IMMOBILE *****************************************************/
/***********************************                   *****************************************************/

// Functors for computing models
struct fn_immob_scale_aqua {
    inline double operator() (double por_imm) {
        return por_imm;
    }
};

struct fn_immob_scale_sorbed {
    inline double operator() (double por_m, double por_imm, double surface_cond, double rock_density) {
        double phi = por_m/(por_m + por_imm);
    	return (1 - phi) * surface_cond * rock_density;
    }
};

struct fn_immob_surface_cond {
    inline double operator() (double por_m, double por_imm) {
        return 1 - por_m - por_imm;
    }
};


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

void SorptionImmob::init_field_models()
{
    eq_fields_->scale_aqua.set(Model<3, FieldValue<3>::Scalar>::create(fn_immob_scale_aqua(), eq_fields_dual_->immob_porosity_), 0.0);
    eq_fields_->scale_sorbed.set(Model<3, FieldValue<3>::Scalar>::create(
            fn_immob_scale_sorbed(), eq_fields_->porosity, eq_fields_dual_->immob_porosity_, eq_fields_->no_sorbing_surface_cond,
	        eq_fields_->rock_density), 0.0);
    eq_fields_->no_sorbing_surface_cond.set(Model<3, FieldValue<3>::Scalar>::create(
            fn_immob_surface_cond(), eq_fields_->porosity, eq_fields_dual_->immob_porosity_), 0.0);
}
