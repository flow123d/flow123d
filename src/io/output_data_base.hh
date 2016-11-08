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
 * @file    output_data_base.hh
 * @brief   
 */

#ifndef OUTPUT_DATA_BASE_HH_
#define OUTPUT_DATA_BASE_HH_



#include <ostream>
#include <string>
#include <typeinfo>

#include "fields/unit_si.hh"
#include "system/global_defs.h"

using namespace std;

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

    /// Types of VTK value
	typedef enum { VTK_INT8, VTK_UINT8, VTK_INT16, VTK_UINT16, VTK_INT32, VTK_UINT32,
                   VTK_FLOAT32, VTK_FLOAT64
    } VTKValueType;

	/**
	 * Destructor of OutputDataBase
	 */
    virtual ~OutputDataBase() {};

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
    virtual void print_binary_all(ostream &out_stream) = 0;

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

    /// Type of stored data
    VTKValueType vtk_type_;

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

};




#endif /* OUTPUT_DATA_BASE_HH_ */
