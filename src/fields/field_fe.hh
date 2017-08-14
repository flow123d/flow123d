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
 * @file    field_fe.hh
 * @brief   
 */

#ifndef FIELD_FE_HH_
#define FIELD_FE_HH_

#include "petscmat.h"
#include "system/system.hh"
#include "fields/field_algo_base.hh"
#include "fields/vec_seq_double.hh"
#include "fields/fe_value_handler.hh"
#include "mesh/mesh.h"
#include "mesh/point.hh"
#include "mesh/bih_tree.hh"
#include "io/element_data_cache.hh"
#include "fem/dofhandler.hh"
#include "input/factory.hh"

#include <memory>



/**
 * Class representing fields given by finite element approximation.
 *
 */
template <int spacedim, class Value>
class FieldFE : public FieldAlgorithmBase<spacedim, Value>
{
public:
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;
    typedef FieldAlgorithmBase<spacedim, Value> FactoryBaseType;

    /**
     * Default constructor, optionally we need number of components @p n_comp in the case of Vector valued fields.
     */
    FieldFE(unsigned int n_comp=0);

    /**
     * Return Record for initialization of FieldFE that is derived from Abstract
     */
    static const Input::Type::Record & get_input_type();

    /**
     * Setter for the finite element data. The mappings are required for computation of local coordinates.
     * @param dh   Dof handler.
     * @param map1 1D mapping.
     * @param map2 2D mapping.
     * @param map3 3D mapping.
     * @param data Vector of dof values.
     */
    void set_fe_data(std::shared_ptr<DOFHandlerMultiDim> dh,
    		MappingP1<1,3> *map1,
    		MappingP1<2,3> *map2,
    		MappingP1<3,3> *map3,
			VectorSeqDouble *data);

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);

	/**
	 * Initialization from the input interface.
	 */
	virtual void init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data);


    /**
     * Update time and possibly update data from GMSH file.
     */
    bool set_time(const TimeStep &time) override;


    /**
     * Set target mesh.
     */
    void set_mesh(const Mesh *mesh, bool boundary_domain) override;


    /**
     * Copy data vector to given output ElementDataCache
     */
    void fill_data_to_cache(ElementDataCache<double> &output_data_cache);


    /**
     * Return size of vector of data stored in Field
     */
    inline unsigned int data_size() const {
    	return data_vec_->size();
    }


    /// Destructor.
	virtual ~FieldFE();

private:
	/// Create DofHandler object
	void make_dof_handler(const Mesh *mesh);

	/// Interpolate data over all elements of target mesh.
	void interpolate(ElementDataCache<double>::ComponentDataPtr data_vec);

	/// DOF handler object
    std::shared_ptr<DOFHandlerMultiDim> dh_;
    /// Store data of Field
    VectorSeqDouble *data_vec_;
    /// Array of indexes to data_vec_, used for get/set values
    std::vector<unsigned int> dof_indices;

    /// Value handler that allows get value of 1D elements.
    FEValueHandler<1, spacedim, Value> value_handler1_;
    /// Value handler that allows get value of 2D elements.
    FEValueHandler<2, spacedim, Value> value_handler2_;
    /// Value handler that allows get value of 3D elements.
    FEValueHandler<3, spacedim, Value> value_handler3_;

    /**
     * Used in DOFHandler::distribute_dofs method. Represents 1D element.
     *
     * For correct functionality must be created proper descendant of FiniteElement class.
     */
    FiniteElement<1,3> *fe1_;
    /// Same as previous, but represents 2D element.
    FiniteElement<2,3> *fe2_;
    /// Same as previous, but represents 3D element.
    FiniteElement<3,3> *fe3_;

	/// mesh reader file
	FilePath reader_file_;

	/// field name read from input
	std::string field_name_;

	/// Field flags.
	FieldFlag::Flags flags_;

    /// Registrar of class to factory
    static const int registrar;
};


#endif /* FIELD_FE_HH_ */
