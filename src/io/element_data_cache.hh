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

#include <boost/exception/info.hpp>           // for error_info::~error_info...
#include <memory>                             // for shared_ptr
#include <sstream>                            // for basic_ostream::write
#include <string>                             // for string, operator<<
#include <vector>                             // for vector
#include <armadillo>
#include "io/element_data_cache_base.hh"      // for ElementDataCacheBase
#include "system/exceptions.hh"               // for ExcStream, operator<<
#include "system/index_types.hh"              // for LongIdx
class Tokenizer;
class Distribution;
struct MeshDataHeader;


/// Return type of method that checked data stored in ElementDataCache (NaN values, limits)
typedef enum  {
	ok,              ///< All values are not NaN and are in limits.
	out_of_limits,   ///< Some value(s) is out of limits
	not_a_number     ///< Some value(s) is set to NaN
} CheckResult;


template <typename T>
class ElementDataCache : public ElementDataCacheBase {
/**
 * Container of the field data on elements used as a common data storage for
 * output of various fields using various output formats and to cache data of several fields when reading the input file.
 * This container also perform serialization for the serial output.
 * Read of values from tokenizer and output of values to stream is implemented as it depends on the value type T.
 */
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
     * @param field_name    Field name thas is read
     * @param time          Actual time of data
	 * @param size_of_cache Count of columns of data cache
	 * @param row_vec_size  Count of rows of data cache
	 */
	ElementDataCache(std::string field_name, double time, unsigned int size_of_cache, unsigned int row_vec_size);

    /**
     * \brief Constructor of output ElementDataCache (allow write data)
     *
     * Has fix size of cache.
     *
     * @param field_name Field name is written as parameter to output stream
     * @param n_comp     Given from shape of field
     * @param size       Count of rows of data cache
     */
	ElementDataCache(std::string field_name, unsigned int n_comp, unsigned int size);

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

    void print_yaml_subarray(ostream &out_stream, unsigned int precision, unsigned int begin, unsigned int end) override;

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
     * Make full check of data stored in cache.
     *
     * Method iterates through data and
     *  - checks NaN data values, default_val replaces NaN
     *  - if default_val==NaN and some value(s) is not replaced with valid value return CheckResult::nan
     *  - if some value(s) is out of limits )lower_bound, upper_bound) return CheckResult::out_of_limits
     *  - in other cases return CheckResult::ok
     *
     * Method is executed only once.
     */
    CheckResult check_values(double default_val, double lower_bound, double upper_bound);

    /**
     * Scale data vector of given 'component_idx' with scale 'coef'.
     *
     * Method is executed only once and should be called after check_values method.
     * Method can be used e. g. for convert between units.
     */
    void scale_data(double coef);

    /// Implements ElementDataCacheBase::gather.
    std::shared_ptr< ElementDataCacheBase > gather(Distribution *distr, LongIdx *local_to_global) override;

    /// Implements ElementDataCacheBase::element_node_cache_fixed_size.
    std::shared_ptr< ElementDataCacheBase > element_node_cache_fixed_size(std::vector<unsigned int> &offset_vec) override;

    /// Implements ElementDataCacheBase::element_node_cache_optimize_size.
    std::shared_ptr< ElementDataCacheBase > element_node_cache_optimize_size(std::vector<unsigned int> &offset_vec) override;

    /// Implements ElementDataCacheBase::compute_node_data.
    std::shared_ptr< ElementDataCacheBase > compute_node_data(std::vector<unsigned int> &conn_vec, unsigned int data_size) override;

    /// Access i-th element in the data vector of 0th component.
    T& operator[](unsigned int i);

protected:
    /// Allow to hold sign, if data in cache is checked and scale (both can be executed only once)
	enum CheckScaleData {
	    none,      ///< Data is neither checked nor scaled.
		check,     ///< Data is only checked.
	    scale      ///< Data is scaled.
	};


	/// Return MPI data type corresponding with template parameter of cache. Needs template specialization.
    MPI_Datatype mpi_data_type();

	/// Sign, if data in cache is checked and scale.
	CheckScaleData check_scale_data_;

	/**
	 * Table of element data.
	 *
	 * For every components contains vector of element data.
	 */
	CacheData data_;
    
};


#endif /* ELEMENT_DATA_CACHE_HH_ */
