/*
 * output_data.cc
 *
 *  Created on: Jun 28, 2016
 *      Author: jb
 */

#include "io/output_data.hh"
#include "io/output_time.hh"
#include "system/armadillo_tools.hh"

template <typename T>
OutputData<T>::OutputData(std::string field_name, unsigned int n_rows, unsigned int n_cols, unsigned int size)
{
	this->set_vtk_type<T>();
    this->field_name = field_name;
    this->output_field_name = this->field_name;

    this->n_values = size;

    if (n_cols == 1) {
        if (n_rows == 1) {
            this->n_elem_ = N_SCALAR;
            this->n_rows = 1;
            this->n_cols = 1;
        } else {
            if (n_rows > 1) {
                if (n_rows > 3) {
                    xprintf(PrgErr,
                            "Do not support output of vectors with fixed size >3. Field: %s\n",
                            this->field_name.c_str());
                } else {
                    this->n_elem_ = N_VECTOR;
                    this->n_rows = 3;
                    this->n_cols = 1;
                }
            } else {
                THROW(OutputTime::ExcOutputVariableVector() << OutputTime::EI_FieldName(this->field_name));
            }
        }
    } else {
        this->n_elem_ = N_TENSOR;
        this->n_rows = 3;
        this->n_cols = 3;
    }

    this->data_ = ElementDataCache<T>::create_data_cache(1, this->n_values * this->n_elem_);
    arr_ptr = (T *)malloc( this->n_rows * this->n_cols * sizeof(T));
}

/**
 * \brief Destructor of OutputData
 */
template <typename T>
OutputData<T>::~OutputData()
{
	free(arr_ptr);
}


/**
 * Output data element on given index @p idx. Method for writing data
 * to output stream.
 *
 * \note This method is used only by MSH file format.
 */
template <typename T>
void OutputData<T>::print_ascii(ostream &out_stream, unsigned int idx)
{
	ASSERT_LT(idx, this->n_values).error();
	std::vector<T> &vec = *( this->data_[0].get() );
	for(unsigned int i = n_elem_*idx; i < n_elem_*(idx+1); ++i )
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
void OutputData<T>::print_ascii_all(ostream &out_stream)
{
    std::vector<T> &vec = *( this->data_[0].get() );
	for(unsigned int idx = 0; idx < this->n_values; idx++) {
    	for(unsigned int i = n_elem_*idx; i < n_elem_*(idx+1); ++i )
    		out_stream << vec[i] << " ";
    }
}


/// Prints the whole data vector into stream.
template <typename T>
void OutputData<T>::print_binary_all(ostream &out_stream, bool print_data_size)
{
	if (print_data_size) {
		// write size of data
		unsigned long long int data_byte_size = this->n_values * n_elem_ * sizeof(T);
		out_stream.write(reinterpret_cast<const char*>(&data_byte_size), sizeof(unsigned long long int));
	}
    // write data
	std::vector<T> &vec = *( this->data_[0].get() );
    for(unsigned int idx = 0; idx < this->n_values; idx++) {
    	for(unsigned int i = n_elem_*idx; i < n_elem_*(idx+1); ++i )
    		out_stream.write(reinterpret_cast<const char*>(&(vec[i])), sizeof(T));
    }
}


template <typename T>
void OutputData<T>::print_all_yaml(ostream &out_stream, unsigned int precision)
{
    out_stream << "[ ";
	std::vector<T> &vec = *( this->data_[0].get() );
    for(unsigned int idx = 0; idx < this->n_values; idx++) {
        if (idx != 0) out_stream << ", ";
        //typename Value::return_type value;
        //out_stream << field_value_to_yaml( Value::from_raw(value, &vec[n_elem_ * idx]), precision );
        out_stream << field_value_to_yaml( vec[n_elem_ * idx], precision );
    }
    out_stream << " ]";
}


template <typename T>
void OutputData<T>::get_min_max_range(double &min, double &max)
{
	min = std::numeric_limits<double>::max();
	max = std::numeric_limits<double>::min();
	std::vector<T> &vec = *( this->data_[0].get() );
    for(unsigned int idx = 0; idx < this->n_values; idx++) {
    	for(unsigned int i = n_elem_*idx; i < n_elem_*(idx+1); ++i ) {
    		if (vec[i] < min) min = vec[i];
    		if (vec[i] > max) max = vec[i];
    	}
    }
}


/**
 * Store data element of given data value under given index.
 */
template <typename T>
void OutputData<T>::store_value(unsigned int idx, const T * value) {
    operate(idx, value,  [](T& raw, T val) {raw = val;});
};

/**
 * Add value to given index
 */
template <typename T>
void OutputData<T>::add(unsigned int idx, const T * value) {
    operate(idx, value,   [](T& raw, T val) {raw += val;});
};

/**
 * Reset values at given index
 */
template <typename T>
void OutputData<T>::zero(unsigned int idx) {
    operate(idx, arr_ptr,   [](T& raw, T val) {raw = 0;});
};

/**
 * Normalize values at given index
 */
template <typename T>
void OutputData<T>::normalize(unsigned int idx, unsigned int divisor) {
    operate(idx, arr_ptr,   [divisor](T& raw, T val) {raw /= divisor;});
};



// Instantiation of OutputData template.
template class OutputData< double >;
template class OutputData< int >;
template class OutputData< unsigned int >;

