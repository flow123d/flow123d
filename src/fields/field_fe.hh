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
#include "fields/fe_value_handler.hh"
#include "la/vector_mpi.hh"
#include "mesh/mesh.h"
#include "mesh/point.hh"
#include "mesh/bih_tree.hh"
#include "mesh/long_idx.hh"
#include "mesh/range_wrapper.hh"
#include "io/element_data_cache.hh"
#include "io/msh_basereader.hh"
#include "fem/fe_p.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_system.hh"
#include "fem/dofhandler.hh"
#include "fem/finite_element.hh"
#include "fem/dh_cell_accessor.hh"
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
	 * Possible interpolations of input data.
	 */
	enum DataInterpolation
	{
		identic_msh,    //!< identical mesh
		equivalent_msh, //!< equivalent mesh (default value)
		gauss_p0,       //!< P0 interpolation (with the use of Gaussian distribution)
		interp_p0       //!< P0 interpolation (with the use of calculation of intersections)
	};

    /**
     * Default constructor, optionally we need number of components @p n_comp in the case of Vector valued fields.
     */
    FieldFE(unsigned int n_comp=0);

    /**
     * Return Record for initialization of FieldFE that is derived from Abstract
     */
    static const Input::Type::Record & get_input_type();

    /**
     * Return Input selection for discretization type (determines the section of VTK file).
     */
    static const Input::Type::Selection & get_disc_selection_input_type();

    /**
     * Return Input selection that allow to set interpolation of input data.
     */
    static const Input::Type::Selection & get_interp_selection_input_type();

    /**
     * Setter for the finite element data.
     * @param dh              Dof handler.
     * @param data            Vector of dof values, optional (create own vector according to dofhandler).
     * @param component_index Index of component (for vector_view/tensor_view)
     * @return                Data vector of dof values.
     */
    VectorMPI set_fe_data(std::shared_ptr<DOFHandlerMultiDim> dh,
    		unsigned int component_index = 0, VectorMPI dof_values = VectorMPI::sequential(0));

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
    void native_data_to_cache(ElementDataCache<double> &output_data_cache);


    /**
     * Return size of vector of data stored in Field
     */
    unsigned int data_size() const;

    inline std::shared_ptr<DOFHandlerMultiDim> get_dofhandler() const {
    	return dh_;
    }

    inline VectorMPI get_data_vec() const {
    	return data_vec_;
    }

    /// Call begin scatter functions (local to ghost) on data vector
    void local_to_ghost_data_scatter_begin();

    /// Call end scatter functions (local to ghost) on data vector
    void local_to_ghost_data_scatter_end();

    /// Destructor.
	virtual ~FieldFE();

private:
	/// Create DofHandler object
	void make_dof_handler(const Mesh *mesh);

	/// Interpolate data (use Gaussian distribution) over all elements of target mesh.
	void interpolate_gauss(ElementDataCache<double>::ComponentDataPtr data_vec);

	/// Interpolate data (use intersection library) over all elements of target mesh.
	void interpolate_intersection(ElementDataCache<double>::ComponentDataPtr data_vec);

	/// Calculate native data over all elements of target mesh.
	void calculate_native_values(ElementDataCache<double>::ComponentDataPtr data_cache);

	/// Calculate elementwise data over all elements of target mesh.
	void calculate_elementwise_values(ElementDataCache<double>::ComponentDataPtr data_cache);

	/**
	 * Fill data to boundary_dofs_ vector.
	 *
	 * TODO: Temporary solution. Fix problem with merge new DOF handler and boundary Mesh. Will be removed in future.
	 */
	void fill_boundary_dofs();


	/// DOF handler object
    std::shared_ptr<DOFHandlerMultiDim> dh_;
    /// Store data of Field
    VectorMPI data_vec_;
    /// Array of indexes to data_vec_, used for get/set values
    std::vector<LongIdx> dof_indices_;

    /// Value handler that allows get value of 0D elements.
    FEValueHandler<0, spacedim, Value> value_handler0_;
    /// Value handler that allows get value of 1D elements.
    FEValueHandler<1, spacedim, Value> value_handler1_;
    /// Value handler that allows get value of 2D elements.
    FEValueHandler<2, spacedim, Value> value_handler2_;
    /// Value handler that allows get value of 3D elements.
    FEValueHandler<3, spacedim, Value> value_handler3_;

    /**
     * Used in DOFHandler::distribute_dofs method. Represents 0D element.
     *
     * For correct functionality must be created proper descendant of FiniteElement class.
     */
    MixedPtr<FiniteElement> fe_;

	/// mesh reader file
	FilePath reader_file_;

	/// field name read from input
	std::string field_name_;

	/// Specify section where to find the field data in input mesh file
	OutputTime::DiscreteSpace discretization_;

	/// Specify type of FE data interpolation
	DataInterpolation interpolation_;

	/// Field flags.
	FieldFlag::Flags flags_;

    /// Default value of element if not set in mesh data file
    double default_value_;

    /// Accessor to Input::Record
    Input::Record in_rec_;

    /// Is set in set_mesh method. Value true means, that we accept only boundary element accessors in the @p value method.
    bool boundary_domain_;

    /**
     * Hold dofs of boundary elements.
     *
     * TODO: Temporary solution. Fix problem with merge new DOF handler and boundary Mesh. Will be removed in future.
     */
    std::shared_ptr< std::vector<LongIdx> > boundary_dofs_;

    /// Registrar of class to factory
    static const int registrar;
};



/**
 * Method creates FieldFE of existing VectorMPI that represents elementwise field.
 *
 * It's necessary to create new VectorMPI of FieldFE, because DOF handler has generally
 * a different ordering than mesh.
 * Then is need to call fill_output_data method.
 *
 * Temporary solution that will be remove after solving issue 995.
 * TODO: remove duplicate code with 'make_dof_handler'
 */
template <int spacedim, class Value>
std::shared_ptr<FieldFE<spacedim, Value> > create_field(VectorMPI & vec_seq, Mesh & mesh, unsigned int n_comp)
{
	std::shared_ptr<DOFHandlerMultiDim> dh; // DOF handler object allow create FieldFE
	MixedPtr<FiniteElement> fe; // Finite element objects (allow to create DOF handler)

	switch (n_comp) { // prepare FEM objects for DOF handler by number of components
        case 1: { // scalar
            fe = MixedPtr<FE_P_disc>(0);
            break;
        }
        case 3: { // vector
             MixedPtr<FE_P_disc>   fe_base(0) ;
            fe = mixed_fe_system(fe_base, FEType::FEVector, 3);
            break;
        }
        case 9: { // tensor
            MixedPtr<FE_P_disc>   fe_base(0) ;
            fe = mixed_fe_system(fe_base, FEType::FETensor, 9);
            break;
        }
        default:
            ASSERT(false).error("Should not happen!\n");
	}
	// Prepare DOF handler
	DOFHandlerMultiDim dh_par(mesh);
	std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( &mesh, fe);
	dh_par.distribute_dofs(ds);
    dh = dh_par.sequential();

	// Construct FieldFE
	std::shared_ptr< FieldFE<spacedim, Value> > field_ptr = std::make_shared< FieldFE<spacedim, Value> >();
	field_ptr->set_fe_data(dh, 0, VectorMPI::sequential(vec_seq.size()) );
	return field_ptr;
}


/**
 * Fill data to VecSeqDouble in order corresponding with element DOFs.
 *
 * Set data to data vector of field in correct order according to values of DOF handler indices.
 *
 * Temporary solution that will be remove after solving issue 995.
 */
template <int spacedim, class Value>
void fill_output_data(VectorMPI & vec_seq, std::shared_ptr<FieldFE<spacedim, Value> > field_ptr)
{
	auto dh = field_ptr->get_dofhandler();
	unsigned int ndofs = dh->max_elem_dofs();
	unsigned int idof; // iterate over indices
	std::vector<LongIdx> indices(ndofs);

	/*for (auto cell : dh->own_range()) {
		cell.get_loc_dof_indices(indices);
		for(idof=0; idof<ndofs; idof++) field_ptr->get_data_vec()[ indices[idof] ] = (*vec_seq.data_ptr())[ ndofs*cell.elm_idx()+idof ];
	}*/

	// Fill DOF handler of FieldFE with correct permutation of data corresponding with DOFs.
	for (auto ele : dh->mesh()->elements_range()) {
		dh->cell_accessor_from_element(ele.idx()).get_loc_dof_indices(indices);
		for(idof=0; idof<ndofs; idof++) field_ptr->get_data_vec()[ indices[idof] ] = (*vec_seq.data_ptr())[ ndofs*ele.idx()+idof ];
	}
}


#endif /* FIELD_FE_HH_ */
