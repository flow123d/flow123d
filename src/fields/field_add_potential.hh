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
 * @file    field_add_potential.hh
 * @brief   
 */

#ifndef FIELD_ADD_POTENTIAL_HH_
#define FIELD_ADD_POTENTIAL_HH_

#include <armadillo>
#include <memory>

#include "fields/field.hh"
#include "fields/field_algo_base.hh"
#include "fields/field_model.hh"
#include "fields/field_coords.hh"


/*******************************************************************************
 * Functors of FieldModels
 */
using Sclr = double;
using Vect = arma::vec3;

// Functor computing piezo_head_p0
struct fn_add_potential {
	inline Sclr operator() (Vect gravity, Vect coords, Sclr pressure) {
        return arma::dot(gravity, coords) + pressure;
    }
};


/**
 * Factory class (descendant of @p Field<...>::FactoryBase) that is necessary
 * for setting pressure values are piezometric head values.
 */
template <int spacedim, class Value>  // <3, FieldValue<3>::Scalar>
class AddPotentialFactory : public Field<spacedim, Value>::FactoryBase {
public:
    /// Constructor.
    AddPotentialFactory( Field<3, FieldValue<3>::VectorFixed > &gravity, FieldCoords &coords, Field<3, FieldValue<3>::Scalar> &inner_field)
    : gravity_(gravity),
	  coords_(coords),
	  inner_field_(inner_field),
      field_name_(inner_field.input_name())
    {}

    typename Field<spacedim,Value>::FieldBasePtr create_field(Input::Record rec, const FieldCommon &) override {
        Input::AbstractRecord field_a_rec;
        if (rec.opt_val(field_name_, field_a_rec)) {

            FieldAlgoBaseInitData init_data(field_name_, Value::NRows_, UnitSI::dimensionless());
            auto inner_field_ptr = FieldAlgorithmBase<spacedim, Value>::function_factory( field_a_rec, init_data );

            // get domain specification
            Input::Array domain_name_array;
            unsigned int id;
            const RegionDB &region_db = inner_field_.mesh()->region_db();
		    if (rec.opt_val("region", domain_name_array)) {
			    std::vector<string> domain_names = region_db.get_and_check_operands(domain_name_array);
	            inner_field_.set(inner_field_ptr, 0.0, domain_names);

            } else if (rec.opt_val("rid", id)) {
                Region region;
                try {
                    region = region_db.find_id(id);
                } catch (RegionDB::ExcUniqueRegionId &e) {
                    e << rec.ei_address();
                    throw;
                }
                if (region.is_valid())
                    inner_field_.set(inner_field_ptr, 0.0, { region.label() });
                else
                    THROW(RegionDB::ExcUnknownRegion() << RegionDB::EI_ID(id) );
            } else {
                inner_field_.set(inner_field_ptr, 0.0); // set on all regions
            }

           	return Model<3, FieldValue<3>::Scalar>::create(fn_add_potential(), gravity_, coords_, inner_field_);

        } else {
            return typename Field<spacedim,Value>::FieldBasePtr();
        }
    }

    bool is_active_field_descriptor(const Input::Record &in_rec, FMT_UNUSED const std::string &input_name) override {
        return in_rec.find<Input::AbstractRecord>(field_name_);
    }

    Field<3, FieldValue<3>::VectorFixed > &gravity_;
    FieldCoords &coords_;
    Field<3, FieldValue<3>::Scalar> &inner_field_;
    std::string field_name_;
};


/**
 * This field is meant to be used to implement two possibilities for initialization of pressure fields in
 * Darcy flows. You can either use pressure of piezo-metric  head which are  related by adding gravity potential.
 * For various reasons we use pressure as the primary variable, so if the user enters piezo-head we need to add the potential to
 * the field he/she has provided. This is done by this class. Unfortunately it introduce one more level of indirection,
 * namely one more virtual call for getting the field value.
 *
 * - The field can not be initialized form the input.
 * - We allow only Scalar Value with element_type double.
 */
template <int spacedim, class Value>
class FieldAddPotential : public FieldAlgorithmBase<spacedim, Value> {
public:
	typedef typename Field<spacedim, Value>::FactoryBase FactoryBaseType;
    typedef typename Space<spacedim>::Point Point;
    /**
     *
     */
    FieldAddPotential( const arma::vec::fixed<spacedim+1> &potential_grad, const Input::AbstractRecord &rec, unsigned int n_comp=0);

    /**
     * Constructor allows to set existing Field as inner_field_
     */
    FieldAddPotential( const arma::vec::fixed<spacedim+1> &potential_grad, std::shared_ptr< FieldAlgorithmBase<spacedim, Value> > inner_field,
    		unsigned int n_comp=0);


    /**
     * Factory class (descendant of @p Field<...>::FactoryBase) that is necessary
     * for setting pressure values are piezometric head values.
     */
    class FieldFactory : public FactoryBaseType {
    public:
    	/// Constructor.
    	FieldFactory(arma::vec::fixed<spacedim+1> potential, std::string field_name)
    	: potential_(potential),
    	  field_name_(field_name)
    	{}

    	typename Field<spacedim,Value>::FieldBasePtr create_field(Input::Record rec, const FieldCommon &) override {
       		Input::AbstractRecord field_a_rec;
        	if (rec.opt_val(field_name_, field_a_rec)) {
        		return std::make_shared< FieldAddPotential<3, FieldValue<3>::Scalar > >( potential_, field_a_rec);
        	} else {
        		return typename Field<spacedim,Value>::FieldBasePtr();
        	}
    	}

    	bool is_active_field_descriptor(const Input::Record &in_rec, FMT_UNUSED const std::string &input_name) override {
    		return in_rec.find<Input::AbstractRecord>(field_name_);
    	}

    	arma::vec::fixed<spacedim+1> potential_;
    	std::string field_name_;
    };


    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);


    /**
     * Update time and possibly update data.
     */
    bool set_time(const TimeStep &time) override;
    
    /// Implements @p FieldAlgirithmBase::set_mesh.
    void set_mesh(const Mesh *mesh, bool boundary_domain) override;
    
    virtual ~FieldAddPotential();

private:
    /// Field to which we add linear potential.
    std::shared_ptr< FieldAlgorithmBase<spacedim, Value> > inner_field_;
    /// Potential gradient.
    arma::vec::fixed<spacedim> grad_;
    /// Potential constant term.
    double zero_level_;
};


#endif /* FIELD_ADD_POTENTIAL_HH_ */
