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
 * @file    element_data_cache.hh
 * @brief   
 */

#ifndef ELEMENT_DATA_CACHE_HH_
#define ELEMENT_DATA_CACHE_HH_

#include <vector>
#include <string>
#include <memory>
#include <ostream>
#include <typeinfo>
#include "system/system.hh"
#include "system/tokenizer.hh"
#include "system/global_defs.h"


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
};



struct MeshDataHeader;


template <typename T>
class ElementDataCache : public ElementDataCacheBase {
public:
	typedef std::shared_ptr< std::vector<T> > ComponentDataPtr;
	typedef std::vector< ComponentDataPtr > CacheData;

	/// Default constructor
	ElementDataCache();

	/**
	 * \brief Constructor of input ElementDataCache (allow read data from GMSH or VTK file)
	 *
	 * Allows set variable size of cache.
	 *
	 * @param data_header   Set data members time_ and field_name_
	 * @param size_of_cache Count of columns of data cache
	 * @param row_vec_size  Count of rows of data cache
	 */
	ElementDataCache(MeshDataHeader data_header, unsigned int size_of_cache, unsigned int row_vec_size);

    /**
     * \brief Constructor of output ElementDataCache (allow write data)
     *
     * Has fix size of cache.
     *
     * @param field_name Field name is written as parameter to output stream
     * @param n_rows     Given from shape of field
     * @param n_cols     Given from shape of field
     * @param size       Count of rows of data cache
     */
	ElementDataCache(std::string field_name, unsigned int n_rows, unsigned int n_cols, unsigned int size);

    /**
     * \brief Destructor of ElementDataCache
     */
    virtual ~ElementDataCache() override;

    /// Return vector of element data for get component.
	ComponentDataPtr get_component_data(unsigned int component_idx);

	/**
	 * Create data cache with given count of columns (\p size_of_cache) and rows (\p row_vec_size).
	 */
	static CacheData create_data_cache(unsigned int size_of_cache, unsigned int row_vec_size);

	/// Implements @p ElementDataCacheBase::read_ascii_data.
	void read_ascii_data(Tokenizer &tok, unsigned int n_components, unsigned int i_row) override;

	/// Implements @p ElementDataCacheBase::read_binary_data.
	void read_binary_data(std::istream &data_stream, unsigned int n_components, unsigned int i_row) override;

    /**
     * Output data element on given index @p idx. Method for writing data
     * to output stream.
     *
     * \note This method is used only by MSH file format.
     */
    void print_ascii(ostream &out_stream, unsigned int idx) override;

    /**
     * \brief Print all data stored in output data ro ascii format
     *
     * TODO: indicate if the tensor data are output in column-first or raw-first order
     *       and possibly implement transposition. Set such property for individual file formats.
     *       Class OutputData stores always in raw-first order.
     */
    void print_ascii_all(ostream &out_stream) override;

    /**
     * \brief Print all data stored in output data to appended binary format
     */
    void print_binary_all(ostream &out_stream, bool print_data_size = true) override;

    void print_all_yaml(ostream &out_stream, unsigned int precision) override;

    /**
     * Store data element of given data value under given index.
     */
    void store_value(unsigned int idx, const T * value);

    /**
     * Add value to given index
     */
    void add(unsigned int idx, const T * value);

    /**
     * Reset values at given index
     */
    void zero(unsigned int idx);

    /**
     * Normalize values at given index
     */
    void normalize(unsigned int idx, unsigned int divisor);

    /**
     * Find minimal and maximal range of stored data
     */
    void get_min_max_range(double &min, double &max) override;

    /**
     * Declaration of new exception info used in following exception
     */
    TYPEDEF_ERR_INFO(EI_FieldName, std::string);

    /**
     * Declaration of exception
     */
    DECLARE_EXCEPTION(ExcOutputVariableVector, << "Can not output field " << EI_FieldName::qval
            << " returning variable size vectors. Try convert to MultiField.\n");

protected:
	/**
	 * Table of element data.
	 *
	 * For every components contains vector of element data.
	 */
	CacheData data_;
    
};

#endif /* ELEMENT_DATA_CACHE_HH_ */
