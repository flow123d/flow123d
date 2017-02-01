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
#include "mesh/mesh.h"
#include "mesh/point.hh"
#include "mesh/bih_tree.hh"
#include "mesh/ngh/include/ngh_interface.hh"
#include "fem/dofhandler.hh"
#include "fem/mapping.hh"
#include "fem/fe_p.hh"
#include "input/factory.hh"



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
    void set_fe_data(DOFHandlerMultiDim *dh,
    		Mapping<1,3> *map1,
    		Mapping<2,3> *map2,
    		Mapping<3,3> *map3,
    		const Vec *data);

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


    /// Destructor.
	virtual ~FieldFE();

private:

    DOFHandlerMultiDim *dh_;
    double *data_;
    const Vec *data_vec_;
    unsigned int *dof_indices;

    Mapping<1,3> *map1_;
    Mapping<2,3> *map2_;
    Mapping<3,3> *map3_;

	FE_P_disc<0,1,3> fe1;
	FE_P_disc<0,2,3> fe2;
	FE_P_disc<0,3,3> fe3;

	/// mesh, which is interpolated
	Mesh* source_mesh_;

	/// mesh reader file
	FilePath reader_file_;

	/// field name read from input
	std::string field_name_;

	/// tree of mesh elements
	BIHTree* bih_tree_;

	/// stored index to last computed element
	unsigned int computed_elm_idx_ = numeric_limits<unsigned int>::max();

	/// vector stored suspect elements in calculating the intersection
	std::vector<unsigned int> searched_elements_;

	/// 3D (tetrahedron) element, used for computing intersection
	TTetrahedron tetrahedron_;

	/// 2D (triangle) element, used for computing intersection
	TTriangle triangle_;

	/// 1D (abscissa) element, used for computing intersection
	TAbscissa abscissa_;

	/// 0D (point) element, used for computing intersection
	TPoint point_, found_point_;

    /// Registrar of class to factory
    static const int registrar;
};


#endif /* FIELD_FE_HH_ */
