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
#include "system/system.hh"
#include "system/tokenizer.hh"
#include "io/output_data_base.hh"


class ElementDataCacheBase {
public:
	/// Constructor.
	ElementDataCacheBase()
	: time_(-std::numeric_limits<double>::infinity()),
	  quantity_name_("") {}

	/// Destructor
	virtual ~ElementDataCacheBase() {}

	/// Getter for time of cache
	double get_time()
	{ return time_; }

	/// Getter for quantity name of cache
	std::string get_quantity_name()
	{ return quantity_name_; }

	/// Check if cache stored actual data
	bool is_actual(double time, std::string quantity_name) {
		return (time_ == time) && (quantity_name_ == quantity_name);
	}

	/**
	 * Read ascii data of given \p i_row from tokenizer
	 */
	virtual void read_ascii_data(Tokenizer &tok, unsigned int n_components, unsigned int i_row)=0;

	/**
	 * Read binary data of given \p i_row from data stream
	 */
	virtual void read_binary_data(std::istream &data_stream, unsigned int n_components, unsigned int i_row)=0;

protected:
	/// time step stored in cache
	double time_;
	/// name of quantity stored in cache
	std::string quantity_name_;
};



struct MeshDataHeader;


template <typename T>
class ElementDataCache : public ElementDataCacheBase, public OutputDataBase {
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
	 * @param data_header   Set data members time_ and quantity_name_
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
     * Perform given function at given index
     */
    template <class Func>
    void operate(unsigned int idx, const T * val, const Func& func) {
        ASSERT_LT_DBG(idx, this->n_values);
        std::vector<T> &vec = *( this->data_[0].get() );
        unsigned int vec_idx = idx*this->n_elem_;
        for(unsigned int i = 0; i < this->n_elem_; i++) {
        	func(vec[vec_idx], val[i]);
        	vec_idx++;
        }
    };


	/**
	 * Table of element data.
	 *
	 * For every components contains vector of element data.
	 */
	CacheData data_;
    
    /**
     * Helper data member represents field.
     *
     * Used in zero and normalize methods.
     */
    T * arr_ptr;

};

#endif /* ELEMENT_DATA_CACHE_HH_ */
