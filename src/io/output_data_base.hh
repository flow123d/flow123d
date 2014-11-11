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
	enum ValueType {
		scalar=1,
		vector=3,
		tensor=9
	};

    virtual ~OutputDataBase() {};
    virtual void print(ostream &out_stream, unsigned int idx) = 0;


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
    ValueType n_elem_;


};




#endif /* OUTPUT_DATA_BASE_HH_ */
