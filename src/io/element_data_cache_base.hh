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

#include "system/tokenizer.hh"
#include "system/system.hh"

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
	  field_name_("") {}

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
    virtual void print_ascii_all(ostream &out_stream) = 0;

    /**
     * Print all data in binary format at once stored in database
     */
    virtual void print_binary_all(ostream &out_stream, bool print_data_size = true) = 0;

    /**
     * Print stored values in the YAML format (using JSON like arrays).
     * Used for output of observe values.
     */
    virtual void print_all_yaml(ostream &out_stream, unsigned int precision) = 0;

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
    inline std::string field_units() {
    	return this->field_units_;
    }

    /**
     * Get number of data values.
     */
    inline unsigned int n_values() {
    	return this->n_values_;
    }

    /**
     * Get number of data elements per data value.
     */
    inline NumCompValueType n_elem() {
    	return this->n_elem_;
    }

    /// Get type of stored data
    inline VTKValueType vtk_type() {
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
    NumCompValueType n_elem_;

    /// Type of stored data
    VTKValueType vtk_type_;

    /// Hash of DOF handler (attribute of native VTK data)
    std::size_t dof_handler_hash_;
};


#endif /* ELEMENT_DATA_CACHE_BASE_HH_ */
