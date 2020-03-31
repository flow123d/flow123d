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
 * @file    element_data_cache.cc
 * @brief
 */


#include <limits>
#include <ostream>
#include "io/element_data_cache.hh"
#include "io/msh_basereader.hh"
#include "la/distribution.hh"
#include "system/armadillo_tools.hh"
#include "system/system.hh"
#include "system/tokenizer.hh"
#include "boost/lexical_cast.hpp"



template <typename T>
ElementDataCache<T>::ElementDataCache()
: ElementDataCacheBase(),
  check_scale_data_(CheckScaleData::none) {}


template <typename T>
ElementDataCache<T>::ElementDataCache(std::string field_name, double time, unsigned int size_of_cache, unsigned int row_vec_size)
: check_scale_data_(CheckScaleData::none)
{
	this->time_ = time;
	this->field_input_name_ = field_name;
	this->data_ = create_data_cache(size_of_cache, row_vec_size);
}


template <typename T>
ElementDataCache<T>::ElementDataCache(std::string field_name, unsigned int n_comp, unsigned int size)
: check_scale_data_(CheckScaleData::none)
{
	this->set_vtk_type<T>();
    this->field_name_ = field_name;
    this->field_input_name_ = this->field_name_;

    this->n_values_ = size;
    ASSERT_GT(n_comp, 0)(field_name).error("Output field returning variable size vectors. Try convert to MultiField.");
    this->n_comp_ = n_comp;

    this->data_ = ElementDataCache<T>::create_data_cache(1, this->n_values_ * this->n_comp_);
}


template <typename T>
ElementDataCache<T>::~ElementDataCache() {}


template <typename T>
typename ElementDataCache<T>::ComponentDataPtr ElementDataCache<T>::get_component_data(unsigned int component_idx) {
	ASSERT_LT(component_idx, data_.size()).error("Index of component is out of range.\n");
	return data_[component_idx];
}


template <typename T>
typename ElementDataCache<T>::CacheData ElementDataCache<T>::create_data_cache(unsigned int size_of_cache, unsigned int row_vec_size) {
    typename ElementDataCache<T>::CacheData data_cache(size_of_cache);
    for (unsigned int i=0; i<size_of_cache; ++i) {
		typename ElementDataCache<T>::ComponentDataPtr row_vec = std::make_shared<std::vector<T>>();
		row_vec->resize(row_vec_size, numeric_limits<T>::signaling_NaN());
		data_cache[i] = row_vec;
    }

    return data_cache;
}


template <typename T>
void ElementDataCache<T>::read_ascii_data(Tokenizer &tok, unsigned int n_components, unsigned int i_row) {
	unsigned int idx;
	for (unsigned int i_vec=0; i_vec<data_.size(); ++i_vec) {
		idx = i_row * n_components;
		std::vector<T> &vec = *( data_[i_vec].get() );
		for (unsigned int i_col=0; i_col < n_components; ++i_col, ++idx) {
			vec[idx] = boost::lexical_cast<T>(*tok);
			++tok;
		}
	}
}


template <typename T>
void ElementDataCache<T>::read_binary_data(std::istream &data_stream, unsigned int n_components, unsigned int i_row) {
	unsigned int idx;
	for (unsigned int i_vec=0; i_vec<data_.size(); ++i_vec) {
		idx = i_row * n_components;
		std::vector<T> &vec = *( data_[i_vec].get() );
		for (unsigned int i_col=0; i_col < n_components; ++i_col, ++idx) {
			data_stream.read(reinterpret_cast<char *>(&vec[idx]), sizeof(T));
		}
	}
}


/**
 * Output data element on given index @p idx. Method for writing data
 * to output stream.
 *
 * \note This method is used only by MSH file format.
 */
template <typename T>
void ElementDataCache<T>::print_ascii(ostream &out_stream, unsigned int idx)
{
	ASSERT_LT(idx, this->n_values_).error();
	std::vector<T> &vec = *( this->data_[0].get() );
	for(unsigned int i = n_comp_*idx; i < n_comp_*(idx+1); ++i )
		out_stream << vec[i] << " ";
}

/**
 * \brief Print all data stored in output data
 *
 * TODO: indicate if the tensor data are output in column-first or raw-first order
 *       and possibly implement transposition. Set such property for individual file formats.
 *       Class OutputData stores always in raw-first order.
 */
template <typename T>
void ElementDataCache<T>::print_ascii_all(ostream &out_stream)
{
    std::vector<T> &vec = *( this->data_[0].get() );
	for(unsigned int idx = 0; idx < this->n_values_; idx++) {
    	for(unsigned int i = n_comp_*idx; i < n_comp_*(idx+1); ++i )
    		out_stream << vec[i] << " ";
    }
}


/// Prints the whole data vector into stream.
template <typename T>
void ElementDataCache<T>::print_binary_all(ostream &out_stream, bool print_data_size)
{
	if (print_data_size) {
		// write size of data
		unsigned long long int data_byte_size = this->n_values_ * n_comp_ * sizeof(T);
		out_stream.write(reinterpret_cast<const char*>(&data_byte_size), sizeof(unsigned long long int));
	}
    // write data
	std::vector<T> &vec = *( this->data_[0].get() );
    for(unsigned int idx = 0; idx < this->n_values_; idx++) {
    	for(unsigned int i = n_comp_*idx; i < n_comp_*(idx+1); ++i )
    		out_stream.write(reinterpret_cast<const char*>(&(vec[i])), sizeof(T));
    }
}


template <typename T>
void ElementDataCache<T>::print_yaml_subarray(ostream &out_stream, unsigned int precision, unsigned int begin, unsigned int end)
{
    out_stream << "[ ";
	std::vector<T> &vec = *( this->data_[0].get() );
    for(unsigned int idx = begin; idx < end; idx++) {
        if (idx != begin) out_stream << " , ";
        unsigned int vec_pos = n_comp_ * idx; // position of element value in data cache
        switch (this->n_comp_) {
            case NumCompValueType::N_SCALAR: {
                out_stream << field_value_to_yaml( vec[vec_pos], precision );
                break;
            }
            case NumCompValueType::N_VECTOR: {
                typename arma::Col<T>::template fixed<3> vec_val;
                for (unsigned int i=0; i<3; ++i, ++vec_pos)
                    vec_val(i) = vec[vec_pos];
                out_stream << field_value_to_yaml( vec_val, precision );
                break;
            }
            case NumCompValueType::N_TENSOR: {
                typename arma::Mat<T>::template fixed<3,3> mat_val;
                for (unsigned int i=0; i<3; ++i)
                    for (unsigned int j=0; j<3; ++j, ++vec_pos)
                    	mat_val(i,j) = vec[vec_pos];
                out_stream << field_value_to_yaml( mat_val, precision );
                break;
            }
        }
    }
    out_stream << " ]";
}


template <typename T>
void ElementDataCache<T>::get_min_max_range(double &min, double &max)
{
	min = std::numeric_limits<double>::max();
	max = std::numeric_limits<double>::min();
	std::vector<T> &vec = *( this->data_[0].get() );
    for(unsigned int idx = 0; idx < this->n_values_; idx++) {
    	for(unsigned int i = n_comp_*idx; i < n_comp_*(idx+1); ++i ) {
    		if (vec[i] < min) min = vec[i];
    		if (vec[i] > max) max = vec[i];
    	}
    }
}


/**
 * Store data element of given data value under given index.
 */
template <typename T>
void ElementDataCache<T>::store_value(unsigned int idx, const T * value) {
    ASSERT_LT_DBG(idx, this->n_values_)(this->field_name_);
    std::vector<T> &vec = *( this->data_[0].get() );
    unsigned int vec_idx = idx*this->n_comp_;
    for(unsigned int i = 0; i < this->n_comp_; i++, vec_idx++) {
    	vec[vec_idx] = value[i];
    }
}

/**
 * Add value to given index
 */
template <typename T>
void ElementDataCache<T>::add(unsigned int idx, const T * value) {
    ASSERT_LT_DBG(idx, this->n_values_);
    std::vector<T> &vec = *( this->data_[0].get() );
    unsigned int vec_idx = idx*this->n_comp_;
    for(unsigned int i = 0; i < this->n_comp_; i++, vec_idx++) {
    	vec[vec_idx] += value[i];
    }
}

/**
 * Reset values at given index
 */
template <typename T>
void ElementDataCache<T>::zero(unsigned int idx) {
    ASSERT_LT_DBG(idx, this->n_values_);
    std::vector<T> &vec = *( this->data_[0].get() );
    unsigned int vec_idx = idx*this->n_comp_;
    for(unsigned int i = 0; i < this->n_comp_; i++, vec_idx++) {
    	vec[vec_idx] = 0;
    }
}

/**
 * Normalize values at given index
 */
template <typename T>
void ElementDataCache<T>::normalize(unsigned int idx, unsigned int divisor) {
    ASSERT_LT_DBG(idx, this->n_values_);
    std::vector<T> &vec = *( this->data_[0].get() );
    unsigned int vec_idx = idx*this->n_comp_;
    for(unsigned int i = 0; i < this->n_comp_; i++, vec_idx++) {
    	vec[vec_idx] /= divisor;
    }
}

template <typename T>
CheckResult ElementDataCache<T>::check_values(double default_val, double lower_bound, double upper_bound) {
    if (check_scale_data_ != CheckScaleData::none) return CheckResult::ok; // method is executed only once
    check_scale_data_ = CheckScaleData::check;

    bool is_nan = false, out_of_limit = false;
    for (unsigned int j=0; j<data_.size(); ++j) {
        std::vector<T> &vec = *( this->data_[j].get() );
        for(unsigned int i=0; i<vec.size(); ++i) {
            if ( std::isnan(vec[i]) ) {
                if ( std::isnan(default_val) ) is_nan = true;
                else vec[i] = default_val;
            }
            if ( (vec[i] < lower_bound) || (vec[i] > upper_bound) ) out_of_limit = true;
        }
    }

    if (is_nan) return CheckResult::not_a_number;
    else if (out_of_limit) return CheckResult::out_of_limits;
    else return CheckResult::ok;
}

template <typename T>
void ElementDataCache<T>::scale_data(double coef) {
    if (check_scale_data_ == CheckScaleData::scale) return; // method is executed only once
    ASSERT_DBG(check_scale_data_ == CheckScaleData::check).warning("Data should be checked before scaling. Rather call 'check_values'!\n");

    for (unsigned int j=0; j<data_.size(); ++j) {
        std::vector<T> &vec = *( this->data_[j].get() );
        for(unsigned int i=0; i<vec.size(); ++i) {
            vec[i] *= coef;
        }
    }

    check_scale_data_ = CheckScaleData::scale;
}


template <typename T>
std::shared_ptr< ElementDataCacheBase > ElementDataCache<T>::gather(Distribution *distr, LongIdx *local_to_global) {
    std::shared_ptr< ElementDataCache<T> > gather_cache;
    int rank = distr->myp();
    int n_proc = distr->np();

    unsigned int n_global_data;     // global number of data
    int rec_starts[n_proc];         // displacement of first value that is received from each process
    int rec_counts[n_proc];         // number of values that are received from each process
    int *rec_indices_ids = nullptr; // collective values of local to global indexes map of data
    T *rec_data = nullptr;          // collective values of data

    // collects values of data vectors and local to global indexes map on each process
    if (rank==0) {
        for (int i=0; i<n_proc; ++i) {
            rec_starts[i] = distr->begin(i);
            rec_counts[i] = distr->lsize(i);
        }
        n_global_data = distr->size();
        rec_indices_ids = new int [ n_global_data ];
    }
    MPI_Gatherv( local_to_global, distr->lsize(), MPI_INT, rec_indices_ids, rec_counts, rec_starts, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank==0) {
        for (int i=0; i<n_proc; ++i) {
            rec_starts[i] = this->n_comp()*rec_starts[i];
            rec_counts[i] = this->n_comp()*rec_counts[i];
        }
        rec_data = new T [ this->n_comp() * n_global_data ];
    }
    auto &local_cache_vec = *( this->get_component_data(0).get() );
    MPI_Gatherv( &local_cache_vec[0], this->n_comp()*distr->lsize(), this->mpi_data_type(), rec_data, rec_counts, rec_starts, this->mpi_data_type(), 0, MPI_COMM_WORLD) ;

    // create and fill serial cache
    if (rank==0) {
        gather_cache = std::make_shared<ElementDataCache<T>>(this->field_input_name_, (unsigned int)this->n_comp(), n_global_data);
        auto &gather_vec = *( gather_cache->get_component_data(0).get() );
        unsigned int i_global_coord; // counter over serial_mesh->nodes_ cache
        for (unsigned int i=0; i<n_global_data; ++i) {
            i_global_coord = this->n_comp() * rec_indices_ids[i];
            for (unsigned int j=0; j<this->n_comp(); ++j) { //loop over coords
            	ASSERT_LT(i_global_coord+j, gather_vec.size());
                gather_vec[ i_global_coord+j ] = rec_data[ this->n_comp()*i+j ];
            }
        }

        delete[] rec_indices_ids;
        delete[] rec_data;
    }

    return gather_cache;
}


template <typename T>
std::shared_ptr< ElementDataCacheBase > ElementDataCache<T>::element_node_cache_fixed_size(std::vector<unsigned int> &offset_vec) {
    unsigned int n_elem = offset_vec.size();
    std::shared_ptr< ElementDataCache<T> > elem_node_cache = std::make_shared<ElementDataCache<T>>(this->field_input_name_, 4*this->n_comp(), n_elem);
    auto &data_out_vec = *( elem_node_cache->get_component_data(0).get() );
    std::fill( data_out_vec.begin(), data_out_vec.end(), (T)0 );
    auto &data_in_vec = *( this->get_component_data(0).get() );

    unsigned int i_node, i_old, i_new;
    for (unsigned int i_el=0, i_conn=0; i_el<offset_vec.size(); i_el++) {
        for(i_node=4*i_el; i_conn<offset_vec[i_el]; i_conn++, i_node++) {
        	i_old = i_conn*this->n_comp_;
        	i_new = i_node*this->n_comp_;
            for(unsigned int i = 0; i < this->n_comp_; i++) {
            	ASSERT_LT(i_new+i, data_out_vec.size());
            	ASSERT_LT(i_old+i, data_in_vec.size());
            	data_out_vec[i_new+i] = data_in_vec[i_old+i];
            }
        }
    }

    return elem_node_cache;
}


template <typename T>
std::shared_ptr< ElementDataCacheBase > ElementDataCache<T>::element_node_cache_optimize_size(std::vector<unsigned int> &offset_vec) {
    std::shared_ptr< ElementDataCache<T> > elem_node_cache = std::make_shared<ElementDataCache<T>>(this->field_input_name_,
            this->n_comp()/4, offset_vec[offset_vec.size()-1]);
    auto &data_out_vec = *( elem_node_cache->get_component_data(0).get() );
    auto &data_in_vec = *( this->get_component_data(0).get() );

    unsigned int i_node, i_old, i_new;
    for (unsigned int i_el=0, i_conn=0; i_el<offset_vec.size(); i_el++) {
        for(i_node=4*i_el; i_conn<offset_vec[i_el]; i_conn++, i_node++) {
        	i_old = i_node*elem_node_cache->n_comp_;
        	i_new = i_conn*elem_node_cache->n_comp_;
            for(unsigned int i = 0; i < elem_node_cache->n_comp_; i++) {
            	ASSERT_LT(i_new+i, data_out_vec.size());
            	ASSERT_LT(i_old+i, data_in_vec.size());
            	data_out_vec[i_new+i] = data_in_vec[i_old+i];
            }
        }
    }
    return elem_node_cache;
}


template <typename T>
std::shared_ptr< ElementDataCacheBase > ElementDataCache<T>::compute_node_data(std::vector<unsigned int> &conn_vec, unsigned int data_size) {
    ASSERT_EQ(conn_vec.size(), this->n_values());
    unsigned int idx;

    // set output data to zero
    std::shared_ptr< ElementDataCache<T> > node_cache = std::make_shared<ElementDataCache<T>>(this->field_input_name_, this->n_comp(), data_size);
    std::vector<unsigned int> count(data_size, 0);
    for (idx=0; idx < node_cache->n_values(); idx++)
        node_cache->zero(idx);

    auto &data_in_vec = *( this->get_component_data(0).get() );
    for (idx=0; idx < conn_vec.size(); idx++) {
    	ASSERT_LT(conn_vec[idx], node_cache->n_values());
    	ASSERT_LT(this->n_comp()*idx, data_in_vec.size());
    	node_cache->add( conn_vec[idx], &(data_in_vec[this->n_comp() * idx]) );
    	count[ conn_vec[idx] ]++;
    }

    // Compute mean values at nodes
    for(idx=0; idx < node_cache->n_values(); idx++)
    	node_cache->normalize(idx, count[idx]);

    return node_cache;
}


template<>
MPI_Datatype ElementDataCache<double>::mpi_data_type() {
    return MPI_DOUBLE;
}

template<>
MPI_Datatype ElementDataCache<int>::mpi_data_type() {
    return MPI_INT;
}

template<>
MPI_Datatype ElementDataCache<unsigned int>::mpi_data_type() {
    return MPI_UNSIGNED;
}


/// Access i-th element in the data vector.
template <class T>
T& ElementDataCache<T>::operator[](unsigned int i)
{
	std::vector<T> &vec = *( this->data_[0].get() );
    ASSERT_DBG(i < vec.size());
    return vec[i];
}



// explicit instantiation of template class
template class ElementDataCache<unsigned int>;
template class ElementDataCache<int>;
template class ElementDataCache<double>;
