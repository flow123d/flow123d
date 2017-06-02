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
#include "io/element_data_cache.hh"
#include <type_traits>


/**
 * \brief This class is used for storing data that are copied from field.
 *
 *
 */
template <typename T>
class OutputData : public OutputDataBase, public ElementDataCache<T> {
public:

    /**
     * \brief Constructor of templated OutputData
     */
    OutputData(std::string field_name, unsigned int n_rows, unsigned int n_cols, unsigned int size);


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

private:

    /**
     * Perform given function at given index
     */
    template <class Func>
    void operate(unsigned int idx, const T * val, const Func& func) {
        ASSERT_LT_DBG(idx, this->n_values);
        std::vector<T> &vec = *( this->data_[0].get() );
        unsigned int vec_idx = idx*this->n_elem_;
        for(unsigned int i = 0; i < this->n_rows*this->n_cols; i++) {
        	func(vec[vec_idx], val[i]);
        	vec_idx++;
        }
    };


    /**
     * Helper data member represents field.
     *
     * Used in zero and normalize methods.
     */
    T * arr_ptr;


    /**
     * Number of rows and cols in stored data element, valid values are (1,1)
     * for scalar; (3,1) for vectors; (3,3) for tensors
     */
    unsigned int n_rows, n_cols;

};




#endif /* SRC_IO_OUTPUT_DATA_HH_ */
