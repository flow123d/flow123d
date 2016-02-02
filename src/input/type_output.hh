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
 * @file    type_output.hh
 * @brief   
 */

#ifndef TYPE_OUTPUT_HH_
#define TYPE_OUTPUT_HH_


#include "input/input_type_forward.hh"
#include "input/type_record.hh"
#include "input/type_abstract.hh"
#include <string>
#include <ostream>
#include <boost/shared_ptr.hpp>


namespace Input {

namespace Type {

using namespace std;

/// Stores version of program and other base data of application
struct RevNumData {
	std::string version, revision, branch, url;
};

/**
 * @brief Base abstract class for output description of the Input::Type tree.
 *
 * Output into various formats is implemented by derived classes:
 * - OutputText - for human readable description (e. g. used for printout of errors)
 * - OutputJSONMachine - for full machine readable description in standard JSON format
 *
 * Usage example:
 * @code
 * cout << OutputText( &my_record) << endl;
 * @endcode
 *
 * @ingroup input_types
 */
class OutputBase {
public:

    /**
     * @brief Performs output of the documentation into given @p stream. The same effect has the reloaded operator '<<'.
     * Returns reference to the same stream.
     */
    virtual ostream& print(ostream& stream) = 0;

protected:
    /**
     * Types of documentation output.
     *
     * Used in print_impl methods for basic and full printout of Input::Type
     */
    enum DocumentationType {
        key_record,			///< Printout only basic data
        full_record			///< Printout full documentation
    };



    /**
     * Constructor
     *
     * @param type Stores input sequence
     */
    OutputBase();


    /// Destructor
    virtual ~OutputBase();

    // data getters
    /// Gets range of array
    void get_array_sizes(Array array, unsigned int &lower , unsigned int &upper );
    /// Gets description of the given record type.
    const string & get_record_description(const Record *rec);
    /// Gets description of the given abstract type.
    const string & get_abstract_description(const Abstract *a_rec);
    /// Gets range of integer
    void get_integer_bounds(Integer integer, int &lower , int &upper );
    /// Gets range of double
    void get_double_bounds(Double dbl, double &lower , double &upper );
    /// Gets pointer of parent Abstract for given Record
    void get_parent_vec(Record rec, std::vector< boost::shared_ptr<Abstract> > &parent_vec);
    /// Gets pointer of inner type for given Array
    void get_array_type(Array array, boost::shared_ptr<TypeBase> &arr_type);
    /// Gets values of default for given record key
    void get_default(Record::KeyIter it, string &type, string &value);
    /// Gets description of the given selection type.
    const string & get_selection_description(const Selection *sel);
    /// Gets ancestor_.type_name of the given AdHocAbstract type.
    void get_adhoc_parent_name(const AdHocAbstract *a_rec, string &parent_name);
    /// Gets iterator to begin of ancestor_.child_data_ of the given AdHocAbstract type.
    Abstract::ChildDataIter get_adhoc_parent_data(const AdHocAbstract *a_rec);


    /**
     * Perform resolution according to actual @p type (using typeid) and call particular print_impl method.
     */
    void print_base(ostream& stream, const TypeBase *type);


    /**
     * Implements printout of Record @p type
     */
    virtual void print_impl(ostream& stream, const Record *type) = 0;
    /**
     * Implements printout of Array @p type
     */
    virtual void print_impl(ostream& stream, const Array *type) = 0;
    /**
     * Implements printout of Abstract @p type
     */
    virtual void print_impl(ostream& stream, const Abstract *type) = 0;
    /**
     * Implements printout of AdHocAbstract @p type
     */
    virtual void print_impl(ostream& stream, const AdHocAbstract *type) = 0;
    /**
     * Implements printout of Selection @p type
     */
    virtual void print_impl(ostream& stream, const Selection *type) = 0;
    /**
     * Implements printout of Integer @p type
     */
	virtual void print_impl(ostream& stream, const Integer *type) = 0;
    /**
     * Implements printout of Double @p type
     */
	virtual void print_impl(ostream& stream, const Double *type) = 0;
    /**
     * Implements printout of Bool @p type
     */
	virtual void print_impl(ostream& stream, const Bool *type) = 0;
    /**
     * Implements printout of String @p type
     */
	virtual void print_impl(ostream& stream, const String *type) = 0;
    /**
     * Implements printout of FileName @p type
     */
    virtual void print_impl(ostream& stream, const FileName *type) = 0;
    /**
     * Implements printout of FileName @p type
     */
    virtual void print_impl(ostream& stream, const Parameter *type) = 0;

    /**
     * Write out a string with given padding of every new line.
     *
     * @param stream Output stream
     * @param str Printed description
     * @param padding Number of spaces added from left
     * @param hash_count Count of '#' chars in description
     */
    void write_description(std::ostream& stream, const string& str, unsigned int padding, unsigned int hash_count = 1);
    /**
     * Write value stored in @p dft.
     *
     * Enclose value in quotes if it's needed or write info that value is optional or obligatory.
     */
    void write_default_value(std::ostream& stream, Default dft);
	/// Clear all data of processed types
	void clear_processed_types();

    /**
     * Returns true if the type was printed out
     *
     * Checks if the @p processed_types_hash_ contains hash of given type.
     */
    bool was_written(std::size_t hash);


    /// Padding of new level of printout, used where we use indentation.
    static const unsigned int padding_size = 4;
    /// Type of documentation output
    DocumentationType doc_type_;
    /// temporary value for printout of description (used in std::setw function)
    unsigned int size_setw_;

    /// Set of hashes of outputed types. Should replace keys.
    std::set<std::size_t> processed_types_hash_;
    /// Content hash of full IST, value is used for key IST_hash
    TypeBase::TypeHash full_hash_;

};




/**********************************************************************************************************************/

/**
 * @brief Class for create text documentation
 *
 * Record, Abstract and Selection are represented by block of text that contains type name, name, description
 * and count and list of keys (for Record), descendants (for Abstract) or values (for Selection).
 *
 * In list are displayed information about subtypes, e.g. type name, description, value, range of numeric values etc.
 *
 * @ingroup input_types
 */
class OutputText : public OutputBase {
public:
    /**
     * Constructor
     *
     * @param type Stores input sequence
     */
	OutputText(const TypeBase *type) : OutputBase(), type_(type) {}

	ostream& print(ostream& stream) override;
protected:
    void print_impl(ostream& stream, const Record *type);
    void print_impl(ostream& stream, const Array *type);
    void print_impl(ostream& stream, const Abstract *type);
    void print_impl(ostream& stream, const AdHocAbstract *type);
    void print_impl(ostream& stream, const Selection *type);
	void print_impl(ostream& stream, const Integer *type);
	void print_impl(ostream& stream, const Double *type);
	void print_impl(ostream& stream, const Bool *type);
	void print_impl(ostream& stream, const String *type);
    void print_impl(ostream& stream, const FileName *type);
    void print_impl(ostream& stream, const Parameter *type);

    /// Object for which is created printout
    const TypeBase *type_;
};







/**
 * @brief Class for create JSON machine readable documentation
 *
 * Every type is represented by one JSON object, for Selection e.g.:
 *   "name" : (string),
 *   "full_name" : (string),
 *   "type" : "Selection",
 *   "values" : [ { "value" : (int), "name": (string), "description" : (string) } ]
 *
 *   @ingroup input_types
 */
class OutputJSONMachine : public OutputBase {
public:
	OutputJSONMachine(const Record &root_type, RevNumData rev_num_data);

	ostream& print(ostream& stream) override;


protected:
	std::string format_hash(TypeBase::TypeHash hash);
	std::string escape_description(std::string desc);

	/**
	 * Print header of a JSON object describing single TypeBase instance.
	 * @param type - type to print
	 * @param input_type - name of TypeBase descendant.
	 */
	void print_type_header(ostream& stream, const TypeBase *type);

    void print_impl(ostream& stream, const Record *type);
    void print_impl(ostream& stream, const Array *type);
    void print_impl(ostream& stream, const Abstract *type);
    void print_impl(ostream& stream, const AdHocAbstract *type);
    void print_impl(ostream& stream, const Selection *type);
    void print_impl(ostream& stream, const Integer *type);
    void print_impl(ostream& stream, const Double *type);
    void print_impl(ostream& stream, const Bool *type);
    void print_impl(ostream& stream, const String *type);
    void print_impl(ostream& stream, const FileName *type);
    void print_impl(ostream& stream, const Parameter *type);


    /// Print all keys of Abstract type or AdHocAbstract type
    void print_abstract_record_keys(ostream& stream, const Abstract *type);

    /**
     * Print actual version of program.
     */
    void print_program_info(ostream& stream);

    /**
     * Print @p full_hash_ key.
     */
    void print_full_hash(ostream& stream);


    /// Header of the format, printed before call of version print.
    /// see @p print(stream) method
    std::string format_head;
    /// Inner part of the format, printed before first call of recursive print.
    /// see @p print(stream) method
    std::string format_inner;
    /// Tail of the format, printed after all recursive prints are finished and before full hash prints.
    /// see @p print(stream) method
    std::string format_full_hash;
    /// Tail of the format, printed after all recursive prints are finished.
    /// see @p print(stream) method
    std::string format_tail;
    /// Contains version of program and other base data
    RevNumData rev_num_data_;
    
    /// Root type. Have to be printed as first, then we pass through the individual type lists to output everything.
    Record root_type_; 
};


/**
 * Overrides output operator for simple output of the input type tree.
 */
std::ostream& operator<<(std::ostream& stream, OutputText type_output);
std::ostream& operator<<(std::ostream& stream, OutputJSONMachine type_output);



} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_OUTPUT_HH_ */
