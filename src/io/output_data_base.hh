/*
 * output_data_base.hh
 *
 *  Created on: Mar 16, 2014
 *      Author: jb
 */

#ifndef OUTPUT_DATA_BASE_HH_
#define OUTPUT_DATA_BASE_HH_



#include <ostream>
#include <string>

#include "fields/unit_si.hh"

/**
 * \brief Common parent class for templated OutputData.
 *
 * Provides virtual method for output of stored data.
 *
 */
class OutputDataBase {
public:

	/**
	 * Number of components of element data stored in the database.
	 */
	enum NumCompValueType {
		N_SCALAR = 1,
		N_VECTOR = 3,
		N_TENSOR = 9
	};

	/**
	 * Destructor of OutputDataBase
	 */
    virtual ~OutputDataBase() {};

    /**
     * Print one value at given index
     */
    virtual void print(ostream &out_stream, unsigned int idx) = 0;

    /**
     * Print all data at once stored in database
     */
    virtual void print_all(ostream &out_stream) = 0;

    /**
     * Data copied from Field.
     */
    std::string output_field_name;
    std::string field_name;
    UnitSI field_units;

    /**
     * Number of data values.
     */
    unsigned int n_values;

    /**
     * Number of data elements per data value.
     */
    NumCompValueType n_elem_;

};




#endif /* OUTPUT_DATA_BASE_HH_ */
