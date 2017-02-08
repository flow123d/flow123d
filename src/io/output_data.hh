/*
 * output_data.hh
 *
 *  Created on: Jun 28, 2016
 *      Author: jb
 */

#ifndef SRC_IO_OUTPUT_DATA_HH_
#define SRC_IO_OUTPUT_DATA_HH_

#include "io/output_data_base.hh"
#include "system/asserts.hh"
#include <type_traits>


class FieldCommon;

/**
 * \brief This class is used for storing data that are copied from field.
 *
 *
 */
template <class Value>
class OutputData : public OutputDataBase {
public:
    typedef typename Value::element_type ElemType;

    /**
     * \brief Constructor of templated OutputData
     */
    OutputData(const FieldCommon &field, unsigned int size);


    /**
     * \brief Destructor of OutputData
     */
    virtual ~OutputData() override;


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
    void print_binary_all(ostream &out_stream) override;

    void print_all_yaml(ostream &out_stream, unsigned int precision) override;

    /**
     * Store data element of given data value under given index.
     */
    void store_value(unsigned int idx, const Value& value);

    /**
     * Add value to given index
     */
    void add(unsigned int idx, const Value& value);

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

private:

    /**
     * Perform given function at given index
     */
    template <class Func>
    void operate(unsigned int idx, const Value &val, const Func& func) {
        ASSERT_LT_DBG(idx, this->n_values);
        ElemType *ptr = this->data_ + idx*this->n_elem_;
        for(unsigned int i_row = 0; i_row < this->n_rows; i_row++) {
            for(unsigned int i_col = 0; i_col < this->n_cols; i_col++) {
                if (i_row < val.n_rows() && i_col < val.n_cols())
                    func(*ptr, val(i_row, i_col));
                else
                    func(*ptr, 0);
                ptr++;
            }
        }
    };


    /**
     * Computed data values for output stored as continuous buffer of their data elements.
     * One data value has @p n_elem data elements (of type double, int or unsigned int).
     */
    ElemType *data_;


    /**
     * Auxiliary value
     */
    typename Value::return_type aux;


    /**
     * Auxiliary field value envelope over @p aux
     */
    Value val_aux;


    /**
     * Number of rows and cols in stored data element, valid values are (1,1)
     * for scalar; (3,1) for vectors; (3,3) for tensors
     */
    unsigned int n_rows, n_cols;

};




#endif /* SRC_IO_OUTPUT_DATA_HH_ */
