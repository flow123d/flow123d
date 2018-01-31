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
 * @file    reader_internal_csv.hh
 * @brief
 */

#ifndef READER_INTERNAL_CSV_HH_
#define READER_INTERNAL_CSV_HH_


//#include "input/input_type_forward.hh"

//#include "input/storage.hh"
//#include "input/path_base.hh"
#include "input/reader_internal_base.hh"


namespace Input {

using namespace std;


/**
 * @brief Creates storage of part of subtree defined in CSV file.
 *
 * We use this class if input tree contains included CSV file where Array
 * is expected. This type must correspond with type_of_value of Array of Records.
 *
 * @ingroup input
 */
class ReaderInternalCsvInclude : public ReaderInternalBase {
public:
	/// Constructor
	ReaderInternalCsvInclude();

	/**
	 * Create storage of subtree defined in CSV file of given @p array.
	 *
	 * Processing of subtree included CSV file:
	 * 1. We expect Record with special structure and tag 'include_csv' in place of Array.
	 * 2. We load path to CSV file and number of head lines to skip in this file.
	 * 3. We load Record 'format' which determines positions of columns in CSV file.
	 *    Structure of this Record must be equal with subtype of array.
	 * 4. We fill helper @p csv_columns_map_ which allows mapping between columns in CSV
	 *    file and storage of array subtype. We return helper storage of this subtype with
	 *    default values (only determines structure).
	 * 5. We iterate through CSV file. For every line we fill data into helper storage and
	 *    copy this helper storage to storage represents included array.
	 */
    StorageBase * read_storage(PathBase &p, const Type::Array *array);

protected:
    /// List of data types, used for mapping columns in CSV include.
    typedef enum {
    	type_int, type_double, type_bool, type_string, type_sel
    } IncludeDataTypes;

    /// Data of one column of including CSV file.
    struct IncludeCsvData {
    	IncludeDataTypes data_type;
    	vector<unsigned int> storage_indexes;
    	const Type::TypeBase *type;
    };

    StorageBase * make_sub_storage(PathBase &p, const Type::Array *array) override;           ///< Create storage of Type::Array type
    StorageBase * make_sub_storage(PathBase &p, const Type::Selection *selection) override;   ///< Create storage of Type::Selection type
    StorageBase * make_sub_storage(PathBase &p, const Type::Bool *bool_type) override;        ///< Create storage of Type::Bool type
    StorageBase * make_sub_storage(PathBase &p, const Type::Integer *int_type) override;      ///< Create storage of Type::Integer type
    StorageBase * make_sub_storage(PathBase &p, const Type::Double *double_type) override;    ///< Create storage of Type::Double type
    StorageBase * make_sub_storage(PathBase &p, const Type::String *string_type) override;    ///< Create storage of Type::String type

    /// Create vector which contains actual indexes of subtree imported in CSV file.
    vector<unsigned int> create_indexes_vector(PathBase &p);

    /// Set storage of simple input type with value given from CSV file.
    void set_storage_from_csv(unsigned int column_index, StorageBase * item_storage, StorageBase * new_storage);

    /// Checks if value on head represents column position in CSV (starts with '#'). If yes, stores position into \p pos.
    bool check_and_read_position_index(PathBase &p, int &pos);

    /// Depth of CSV included subtree
    unsigned int csv_subtree_depth_;

    /// Map of columns in CSV file to storage of subtree
    map<unsigned int, IncludeCsvData> csv_columns_map_;

};


} // namespace Input

#endif /* READER_INTERNAL_CSV_HH_ */
