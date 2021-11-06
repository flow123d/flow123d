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
 * @file    element_data_cache_base.hh
 * @brief
 */

#ifndef ELEMENT_DATA_CACHE_BASE_HH_
#define ELEMENT_DATA_CACHE_BASE_HH_



#include <ostream>
#include <string>
#include <istream>
#include "system/system.hh"
#include "system/index_types.hh"

class Tokenizer;
class Distribution;


class ElementDataCacheBase {
public:

	/**
	 * Number of components of element data stored in the database.
	 */
	enum NumCompValueType {
		N_SCALAR = 1,
		N_VECTOR = 3,
		N_TENSOR = 9
	};

    /// Types of VTK value
	typedef enum { VTK_INT8, VTK_UINT8, VTK_INT16, VTK_UINT16, VTK_INT32, VTK_UINT32,
                   VTK_FLOAT32, VTK_FLOAT64
    } VTKValueType;

	/// Constructor.
	ElementDataCacheBase()
	: time_(-std::numeric_limits<double>::infinity()),
	  field_name_(""), is_dummy_(false) {}

	/// Destructor
	virtual ~ElementDataCacheBase() {}

	/// Getter for time of cache
	double get_time()
	{ return time_; }

	/// Getter for quantity name of cache
	std::string field_input_name()
	{ return field_input_name_; }

	/// Check if cache stored actual data
	bool is_actual(double time, std::string field_name) {
		return (time_ == time) && (field_input_name_ == field_name);
	}

	/**
	 * Read ascii data of given \p i_row from tokenizer
	 */
	virtual void read_ascii_data(Tokenizer &tok, unsigned int n_components, unsigned int i_row)=0;

	/**
	 * Read binary data of given \p i_row from data stream
	 */
	virtual void read_binary_data(std::istream &data_stream, unsigned int n_components, unsigned int i_row)=0;

    /**
     * Print one value at given index in ascii format
     */
    virtual void print_ascii(ostream &out_stream, unsigned int idx) = 0;

    /**
     * Print all data in ascii format at once stored in database
     */
    virtual void print_ascii_all(ostream &out_stream, unsigned int start=0) = 0;

    /**
     * Print all data in binary format at once stored in database
     */
    virtual void print_binary_all(ostream &out_stream, bool print_data_size = true, unsigned int start = 0) = 0;

    /**
     * Print stored values in the YAML format (using JSON like arrays).
     * Used for output of observe values.
     */
    virtual void print_yaml_subarray(ostream &out_stream, unsigned int precision, unsigned int begin, unsigned int end) = 0;

    /**
     * Find minimal and maximal range of stored data
     */
    virtual void get_min_max_range(double &min, double &max) = 0;

    /**
     * Set string representation of SI units.
     */
    void set_field_units(std::string unit_string) {
    	this->field_units_ = unit_string;
    }

    /**
     * Set string representation of SI units.
     */
    void set_n_values(unsigned int n_values) {
    	this->n_values_ = n_values;
    }

    /**
     * Get string representation of SI units.
     */
    inline std::string field_units() const {
    	return this->field_units_;
    }

    /**
     * Get number of data values.
     */
    inline unsigned int n_values() const {
    	return this->n_values_;
    }

    /**
     * Is only true when the object is DummyElementDataCache.
     */
    inline bool is_dummy() const {
        return this->is_dummy_;
    }

    /**
     * Get number of data elements per data value.
     */
    inline unsigned int n_comp() const {
    	return this->n_comp_;
    }

    /**
     * Get number of DOFs per element.
     */
    inline unsigned int n_dofs_per_element() const {
        return this->n_dofs_per_element_;
    }

    /// Get type of stored data
    inline VTKValueType vtk_type() const {
    	return this->vtk_type_;
    }

    /**
     * Get dof_handler_hash_ value.
     */
    inline std::size_t dof_handler_hash() const {
    	return this->dof_handler_hash_;
    }

    /**
     * Set dof_handler_hash_ value.
     */
    inline void set_dof_handler_hash(std::size_t hash) {
    	this->dof_handler_hash_ = hash;
    }

    /**
     * Get fe_type_ value.
     */
    inline std::string fe_type() const {
    	return this->fe_type_;
    }

    /**
     * Method for gathering parallel data to serial cache.
     *
     * Gather data of individual processes to serial cache that is created only on zero process.
     * Other processes return uninitialized shared pointer.
     *
     * @param distr Collective distribution
     * @param local_to_global Maps local indices to global
     */
    virtual std::shared_ptr< ElementDataCacheBase > gather(Distribution *distr, LongIdx *local_to_global)=0;

    /**
     * Create node data cache of constant data size for each elements.
     *
     * Method must be call on node data cache.
     *
     * Every element is represented of 4 nodes maximal. Method returns cache with size = 4*n_elements*n_components.
     * If dimension of element is less than 3, part of data is not used. This construction of node data cache allow
     * call gather of node data for continuous and discontinuous output meshes.
     *
     * @param offset_vec vector of appropriate offsets (number of nodes) of each elements
     */
    virtual std::shared_ptr< ElementDataCacheBase > element_node_cache_fixed_size(std::vector<unsigned int> &offset_vec)=0;

    /**
     * Inverse method to previous.
     *
     * Must be call on node cache with constant data size for each elements. Return data cache, that corresponds with
     * offset vector.
     *
     * @param offset_vec vector of appropriate offsets (number of nodes) of each elements
     */
    virtual std::shared_ptr< ElementDataCacheBase > element_node_cache_optimize_size(std::vector<unsigned int> &offset_vec)=0;

    /**
     * Create final node data cache.
     *
     * Compute average value of each nodes.
     *
     * @param conn_vec Vector of connectivities, holds node indices to values of this cache
     * @param data_size number of nodes
     */
    virtual std::shared_ptr< ElementDataCacheBase > compute_node_data(std::vector<unsigned int> &conn_vec, unsigned int data_size)=0;

protected:
    template <class T>
    void set_vtk_type() {
    	if ( std::is_same<T, double>::value ) {
    		vtk_type_ = VTK_FLOAT64;
    	} else if ( std::is_same<T, unsigned int>::value ) {
    		vtk_type_ = VTK_UINT32;
    	} else if ( std::is_same<T, int>::value ) {
    		vtk_type_ = VTK_INT32;
    	} else {
    		ASSERT(false).error("Unsupported VTK type");
    	}
    }

	/// time step stored in cache
	double time_;

	/// name of field stored in cache
    std::string field_input_name_;
    std::string field_name_;

    /**
     * Data copied from Field.
     */
    std::string field_units_;

    /**
     * Number of data values.
     */
    unsigned int n_values_;

    /**
     * Number of data elements per data value.
     */
    unsigned int n_comp_;

    /// Type of stored data
    VTKValueType vtk_type_;

    /// Hash of DOF handler (attribute of native VTK data)
    std::size_t dof_handler_hash_;

    /**
     * FiniteElement type (attribute of native VTK data)
     */
    std::string fe_type_;

    /**
     * Number of DOFs per element (attribute of native VTK data)
     */
    unsigned int n_dofs_per_element_;

    /// Is true for DummyElementDataCache
    bool is_dummy_;
};





/**
 * Auxiliary implementation of ElementDataCacheBase that performs output of single zero data for the fields that are
 * off for current time frame.
 */
class DummyElementDataCache : public ElementDataCacheBase {
public:

    DummyElementDataCache(std::string field_name_in, unsigned int n_comp_in)
    {
        this->field_input_name_ = field_name_in;
        this->n_comp_ = n_comp_in;
        this->n_values_ = 1;
        this->is_dummy_ = true;
    }

    virtual ~DummyElementDataCache() override
    {}

    void print_ascii(ostream &out_stream, unsigned int) override
    {
        for(unsigned int i=0; i< n_comp_;i++) out_stream << 0 << " ";
    }

    void print_ascii_all(ostream &out_stream, unsigned int start=0) override
    {
        for(unsigned int i=start; i< n_comp_;i++) out_stream << 0 << " ";
    }

    void print_binary_all(ostream &, bool, unsigned int) override
    {
        ASSERT(false).error("Not implemented.");
    }

    void print_yaml_subarray(ostream &, unsigned int, unsigned int , unsigned int) override
    {}

    void get_min_max_range(double &, double &) override
    {}

    void read_ascii_data(Tokenizer &, unsigned int, unsigned int ) override
    {}

    void read_binary_data(std::istream &, unsigned int, unsigned int) override
    {}

    std::shared_ptr< ElementDataCacheBase > gather(Distribution *, LongIdx *) override
    {
    	return std::make_shared<DummyElementDataCache>(this->field_input_name_, this->n_comp_);
    }

    std::shared_ptr< ElementDataCacheBase > element_node_cache_fixed_size(std::vector<unsigned int> &) override
    {
    	return std::make_shared<DummyElementDataCache>(this->field_input_name_, this->n_comp_);
    }

    std::shared_ptr< ElementDataCacheBase > element_node_cache_optimize_size(std::vector<unsigned int> &) override
    {
    	return std::make_shared<DummyElementDataCache>(this->field_input_name_, this->n_comp_);
    }

    std::shared_ptr< ElementDataCacheBase > compute_node_data(std::vector<unsigned int> &, unsigned int ) override
    {
    	return std::make_shared<DummyElementDataCache>(this->field_input_name_, this->n_comp_);
    }

};

#endif /* ELEMENT_DATA_CACHE_BASE_HH_ */
