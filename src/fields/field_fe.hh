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
#include "system/index_types.hh"
#include "fields/field_algo_base.hh"
#include "fields/fe_value_handler.hh"
#include "la/vector_mpi.hh"
#include "mesh/mesh.h"
#include "mesh/point.hh"
#include "mesh/bih_tree.hh"
#include "mesh/range_wrapper.hh"
#include "io/element_data_cache.hh"
#include "io/msh_basereader.hh"
#include "fem/fe_p.hh"
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
    virtual void value_list (const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);

    /**
     * Overload @p FieldAlgorithmBase::cache_update
     */
    void cache_update(FieldValueCache<typename Value::element_type> &data_cache,
			ElementCacheMap &cache_map, unsigned int region_idx) override;

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

	/// Initialize FEValues object of given dimension.
	template <unsigned int dim>
	Quadrature init_quad(std::shared_ptr<EvalPoints> eval_points);

    Armor::ArmaMat<typename Value::element_type, Value::NRows_, Value::NCols_> handle_fe_shape(unsigned int dim,
            unsigned int i_dof, unsigned int i_qp, unsigned int comp_index);


	/// DOF handler object
    std::shared_ptr<DOFHandlerMultiDim> dh_;
    /// Store data of Field
    VectorMPI data_vec_;

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
    std::shared_ptr< std::vector<IntIdx> > boundary_dofs_;

    /// List of FEValues objects of dimensions 0,1,2,3 used for value calculation
    std::vector<FEValues<spacedim>> fe_values_;

    /// Registrar of class to factory
    static const int registrar;
};


/** Create FieldFE from dhf handler */
template <int spacedim, class Value>
std::shared_ptr<FieldFE<spacedim, Value> > create_field_fe(std::shared_ptr<DOFHandlerMultiDim> dh,
                                                           unsigned int comp = 0,
                                                           VectorMPI *vec = nullptr)
{
	// Construct FieldFE
	std::shared_ptr< FieldFE<spacedim, Value> > field_ptr = std::make_shared< FieldFE<spacedim, Value> >();
    if (vec == nullptr)
	    field_ptr->set_fe_data( dh, comp, dh->create_vector() );
    else
        field_ptr->set_fe_data( dh, comp, *vec );
    
	return field_ptr;
}


/** Create FieldFE with parallel VectorMPI from finite element */
template <int spacedim, class Value>
std::shared_ptr<FieldFE<spacedim, Value> > create_field_fe(Mesh & mesh, const MixedPtr<FiniteElement> &fe)
{
	// Prepare DOF handler
	std::shared_ptr<DOFHandlerMultiDim> dh_par = std::make_shared<DOFHandlerMultiDim>(mesh);
	std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( &mesh, fe);
	dh_par->distribute_dofs(ds);

	return create_field_fe<spacedim,Value>(dh_par);
}


#endif /* FIELD_FE_HH_ */
