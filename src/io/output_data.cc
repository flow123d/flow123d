/*
 * output_data.cc
 *
 *  Created on: Jun 28, 2016
 *      Author: jb
 */

#include "io/output_data.hh"
#include "fields/field_values.hh"
#include "fields/field_common.hh"
#include "io/output_time.hh"
#include "system/armadillo_tools.hh"

template <class Value>
OutputData<Value>::OutputData(const FieldCommon &field, unsigned int size)
: val_aux(aux)
{
	this->set_vtk_type<ElemType>();
    this->field_name = field.name();
    this->field_units = field.units();
    this->output_field_name = this->field_name;

    this->n_values = size;

    if (val_aux.NCols_ == 1) {
        if (val_aux.NRows_ == 1) {
            this->n_elem_ = N_SCALAR;
            this->n_rows = 1;
            this->n_cols = 1;
        } else {
            if (val_aux.NRows_ > 1) {
                if (val_aux.NRows_ > 3) {
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

    data_ = new ElemType[this->n_values * this->n_elem_];
}

/**
 * \brief Destructor of OutputData
 */
template <class Value>
OutputData<Value>::~OutputData()
{
    delete[] this->data_;
}


/**
 * Output data element on given index @p idx. Method for writing data
 * to output stream.
 *
 * \note This method is used only by MSH file format.
 */
template <class Value>
void OutputData<Value>::print_ascii(ostream &out_stream, unsigned int idx)
{
	ASSERT_LT(idx, this->n_values).error();
    ElemType *ptr_begin = this->data_ + n_elem_ * idx;
    for(ElemType *ptr = ptr_begin; ptr < ptr_begin + n_elem_; ptr++ )
        out_stream << *ptr << " ";
}

/**
 * \brief Print all data stored in output data
 *
 * TODO: indicate if the tensor data are output in column-first or raw-first order
 *       and possibly implement transposition. Set such property for individual file formats.
 *       Class OutputData stores always in raw-first order.
 */
template <class Value>
void OutputData<Value>::print_ascii_all(ostream &out_stream)
{
    for(unsigned int idx = 0; idx < this->n_values; idx++) {
        ElemType *ptr_begin = this->data_ + n_elem_ * idx;
        for(ElemType *ptr = ptr_begin; ptr < ptr_begin + n_elem_; ptr++ )
            out_stream << *ptr << " ";
    }
}


/// Prints the whole data vector into stream.
template <class Value>
void OutputData<Value>::print_binary_all(ostream &out_stream)
{
	// write size of data
	unsigned long long int data_byte_size = this->n_values * n_elem_ * sizeof(ElemType);
	out_stream.write(reinterpret_cast<const char*>(&data_byte_size), sizeof(unsigned long long int));
    // write data
    for(unsigned int idx = 0; idx < this->n_values; idx++) {
        ElemType *ptr_begin = this->data_ + n_elem_ * idx;
        for(ElemType *ptr = ptr_begin; ptr < ptr_begin + n_elem_; ptr++ ) {
        	out_stream.write(reinterpret_cast<const char*>(ptr), sizeof(ElemType));
        }
    }
}


template <class Value>
void OutputData<Value>::print_all_yaml(ostream &out_stream, unsigned int precision)
{
    out_stream << "[ ";
    for(unsigned int idx = 0; idx < this->n_values; idx++) {
        if (idx != 0) out_stream << ", ";
        ElemType *ptr_begin = this->data_ + n_elem_ * idx;
        typename Value::return_type value;
        out_stream << field_value_to_yaml( Value::from_raw(value, ptr_begin), precision );
    }
    out_stream << " ]";
}


template <class Value>
void OutputData<Value>::get_min_max_range(double &min, double &max)
{
	min = std::numeric_limits<double>::max();
	max = std::numeric_limits<double>::min();
    for(unsigned int idx = 0; idx < this->n_values; idx++) {
        ElemType *ptr_begin = this->data_ + n_elem_ * idx;
        for(ElemType *ptr = ptr_begin; ptr < ptr_begin + n_elem_; ptr++ ) {
           	if ((double)(*ptr) < min) min = (double)(*ptr);
           	if ((double)(*ptr) > max) max = (double)(*ptr);
        }
    }
}


/**
 * Store data element of given data value under given index.
 */
template <class Value>
void OutputData<Value>::store_value(unsigned int idx, const Value& value) {
    operate(idx, value,  [](ElemType& raw, ElemType val) {raw = val;});
};

/**
 * Add value to given index
 */
template <class Value>
void OutputData<Value>::add(unsigned int idx, const Value& value) {
    operate(idx, value,   [](ElemType& raw, ElemType val) {raw += val;});
};

/**
 * Reset values at given index
 */
template <class Value>
void OutputData<Value>::zero(unsigned int idx) {
    operate(idx, val_aux,   [](ElemType& raw, ElemType val) {raw = 0;});
};

/**
 * Normalize values at given index
 */
template <class Value>
void OutputData<Value>::normalize(unsigned int idx, unsigned int divisor) {
    operate(idx, val_aux,   [divisor](ElemType& raw, ElemType val) {raw /= divisor;});
};



// Instantiation of OutputData template.
template class OutputData< FieldValue<0>::Enum >;
template class OutputData< FieldValue<0>::Integer >;
template class OutputData< FieldValue<0>::Scalar >;

template class OutputData< FieldValue<2>::VectorFixed >;
template class OutputData< FieldValue<2>::TensorFixed >;

template class OutputData< FieldValue<3>::VectorFixed >;
template class OutputData< FieldValue<3>::TensorFixed >;

