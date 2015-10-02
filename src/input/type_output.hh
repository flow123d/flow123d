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


#include "input/type_base.hh"
#include "input/type_record.hh"
#include <boost/regex.hpp>


namespace Input {

namespace Type {

/**
 * @brief Base abstract class for output description of the Input::Type tree.
 *
 * Output into various formats is implemented by derived classes:
 * - OutputText - for human readable description
 * - OutputJSONTemplate - for creating a template of the input file
 * - OutputLatex - for printed documentation (with support of particular Latex macros)
 * - OutputJSONMachine - for full machine readable description in standard JSON format
 *
 * Usage example:
 * @code
 * cout << OutputText( &my_record, 3) << endl;
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
    virtual ostream& print(ostream& stream);

    /**
     * @brief Initializes and allocates regular expression filter @p regex_filter.
     *
     * Full names of Input::Type::Record objects are passed through the filter, deleting the first match of the regular expression given by @p regex_filter.
     * Full Record description is performed only for the first occurrence of the filtered name, further Records with same filtered name are ignored
     * (reported only in short descriptions of individual keys or array subtype, etc.) See @p ProcessedTypes for details.
     */
	void set_filter(string regex_filter);

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
     * @param depth Depth of output
     */
    OutputBase(const TypeBase *type, unsigned int depth = 0);


    /// Destructor
    virtual ~OutputBase();

    // data getters
    /// Gets range of array
    void get_array_sizes(Array array, unsigned int &lower , unsigned int &upper );
    /// Gets description of the given record type.
    const string & get_record_description(const Record *rec);
    /// Gets description of the given abstract type.
    const string & get_abstract_description(const AbstractRecord *a_rec);
    /// Gets record key for given index
    void get_record_key(Record rec, unsigned int key_idx, Record::Key &key);
    /// Gets range of integer
    void get_integer_bounds(Integer integer, int &lower , int &upper );
    /// Gets range of double
    void get_double_bounds(Double dbl, double &lower , double &upper );
    /// Gets pointer of parent AbstractRecord for given Record
    void get_parent_ptr(Record rec, boost::shared_ptr<AbstractRecord> &parent_ptr);
    /// Gets pointer of inner type for given Array
    void get_array_type(Array array, boost::shared_ptr<TypeBase> &arr_type);
    /// Gets values of default for given record key
    void get_default(Record::KeyIter it, string &type, string &value);
    /// Gets description of the given selection type.
    const string & get_selection_description(const Selection *sel);
    /// Gets parent_name_ of the given AdHocAbstractRecord type.
    const string & get_adhoc_parent_name(const AdHocAbstractRecord *a_rec);
    /// Gets iterator to begin of parent_data_ of the given AdHocAbstractRecord type.
    AbstractRecord::ChildDataIter get_adhoc_parent_data(const AdHocAbstractRecord *a_rec);


    /// Gets pointer of inner data for given Record
    const void * get_record_data(const Record *rec);
    /// Gets pointer of inner data for given AbstractRecord
    const void * get_abstract_record_data(const AbstractRecord *a_rec);
    /// Gets pointer of inner data for given Selection
    const void * get_selection_data(const Selection *sel);
    /// Gets pointer of inner data for given Array
    const void * get_array_data(const Array *array);
    /// Gets pointer of inner data for given TypeBase
    const void * get_type_base_data(const TypeBase *type);


    /**
     * Perform resolution according to actual @p type (using typeid) and call particular print_impl method.
     */
    void print(ostream& stream, const TypeBase *type, unsigned int depth);


    /**
     * Implements printout of Record @p type
     */
    virtual void print_impl(ostream& stream, const Record *type, unsigned int depth) = 0;
    /**
     * Implements printout of Array @p type
     */
    virtual void print_impl(ostream& stream, const Array *type, unsigned int depth) = 0;
    /**
     * Implements printout of AbstractRecord @p type
     */
    virtual void print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth) = 0;
    /**
     * Implements printout of AdHocAbstractRecord @p type
     */
    virtual void print_impl(ostream& stream, const AdHocAbstractRecord *type, unsigned int depth) = 0;
    /**
     * Implements printout of Selection @p type
     */
    virtual void print_impl(ostream& stream, const Selection *type, unsigned int depth) = 0;
    /**
     * Implements printout of Integer @p type
     */
	virtual void print_impl(ostream& stream, const Integer *type, unsigned int depth) = 0;
    /**
     * Implements printout of Double @p type
     */
	virtual void print_impl(ostream& stream, const Double *type, unsigned int depth) = 0;
    /**
     * Implements printout of Bool @p type
     */
	virtual void print_impl(ostream& stream, const Bool *type, unsigned int depth) = 0;
    /**
     * Implements printout of String @p type
     */
	virtual void print_impl(ostream& stream, const String *type, unsigned int depth) = 0;
    /**
     * Implements printout of FileName @p type
     */
    virtual void print_impl(ostream& stream, const FileName *type, unsigned int depth) = 0;

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


    /// Padding of new level of printout, used where we use indentation.
    static const unsigned int padding_size = 4;
    /// Object for which is created printout
    const TypeBase *type_;
    /// Depth of printout (for value 0 is printed all input tree)
    unsigned int depth_;
    /// Type of documentation output
    DocumentationType doc_type_;
    /// temporary value for printout of description (used in std::setw function)
    unsigned int size_setw_;

    /// Header of the format, printed before first call of recursive print.
    /// see @p print(stream) method
    std::string format_head;
    /// Tail of the format, printed after all recursive prints are finished.
    /// see @p print(stream) method
    std::string format_tail;

    /**
     * @brief Internal data class.
     * Contains flags of written Input::Types objects and functionality of regular expression filter of Input::Types full names.
     *
     * Flags are stored to struct that contains unique internal data pointer of complex Input::Type,
     * flag if extensive documentation was printed and reference to Input::Type.
     *
     * Regular expression filter is optional and stores printed Input::Type by filtered full_name.
     * Input::Types with similar full names are printed only once.
     */
    class ProcessedTypes {
    public:

        /**
         * Structure for flags about output of one TypeBase object in input tree
         * Stores types what was printed
         */
        struct Key {
        	unsigned int key_index;          	///< Position inside the record.
        	const void * type_data_;            ///< Pointer to internal data of type.
        	mutable string reference_;       	///< Reference to type.
        };
        /**
         * Public typedef of constant iterator into array of keys.
         */
        typedef std::vector<struct Key>::const_iterator KeyIter;


    	/// Clear all data of processed types
    	void clear();

        /**
         * Interface to mapping key -> index. Returns index (in continuous array) for given type.
         */
        unsigned int type_index(const void * type_data) const;

    	/// Destructor, deallocates filter_ if it was allocated.
    	~ProcessedTypes();

        /**
         * Returns true if the type was printed out
         *
         * Checks if the ProcessedTypes contains key of given type and key has true flag extensive_doc_
         * or if the ProcessedTypes contains type of given full_name when regular expression filter_ is initialized.
         */
    	bool was_written(const void * type_data, string full_name);

    	/**
    	 * Marks type as written.
    	 *
    	 * Inserts type to key_to_index map.
    	 * If regular expression filter_ is initialized marks filtered full_name of type as written.
    	 */
    	void mark_written(const void *type_data, string full_name, string reference = "");
        /**
         * Returns reference_ string of key of given type.
         */
        const string get_reference(const void * type_data) const;

        /**
         * Combines was_sritten and mark_written for hashes.
         */
        bool was_written(std::size_t hash)
        {
            bool in_set = ( output_hash.find(hash) != output_hash.end() );
            if (! in_set) output_hash.insert(hash);
            return in_set;
        }


        /// Database of valid keys
        std::map<const void *, unsigned int> key_to_index;
        typedef std::map<const void *, unsigned int>::const_iterator key_to_index_const_iter;

        /// Keys in order as they where declared.
        std::vector<struct Key> keys;

        /// Regex filter for full names.
        boost::regex    *filter_;

        /// Regular expression of filter
        string reg_exp_;

        /// Set of processed types by regular expression and full names
        std::set<string> full_type_names;

        /// Set of hashes of outputed types. Should replace keys.
        std::set<std::size_t> output_hash;
    };

    /// Stores flags and references of processed type
    ProcessedTypes doc_flags_;

};


/**********************************************************************************************************************/

/**
 * @brief Class for create text documentation
 *
 * Record, AbstractRecord and Selection are represented by block of text that contains type name, name, description
 * and count and list of keys (for Record), descendants (for AbstractRecord) or values (for Selection).
 *
 * In list are displayed information about subtypes, e.g. type name, description, value, range of numeric values etc.
 *
 * @ingroup input_types
 */
class OutputText : public OutputBase {
public:
	OutputText(const TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

protected:

    void print_impl(ostream& stream, const Record *type, unsigned int depth);
    void print_impl(ostream& stream, const Array *type, unsigned int depth);
    void print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const AdHocAbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const Selection *type, unsigned int depth);
	void print_impl(ostream& stream, const Integer *type, unsigned int depth);
	void print_impl(ostream& stream, const Double *type, unsigned int depth);
	void print_impl(ostream& stream, const Bool *type, unsigned int depth);
	void print_impl(ostream& stream, const String *type, unsigned int depth);
    void print_impl(ostream& stream, const FileName *type, unsigned int depth);


};







/**
 * @brief Class for create and JSON template documentation
 *
 * Every type is represented by JSON object.
 * Type name, description and other data are displayed as comments.
 * Among other data belongs range of numeric values, size limits of arrays, possible values of selections etc.
 *
 *
 * @ingroup input_types
 */
class OutputJSONTemplate : public OutputBase {
public:
    /**
     * Constructor for output of the input type tree with root @p type.
     * The input type tree is searched by DFS algorithm into @p depth.
     */
	OutputJSONTemplate(TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

	/**
	 * Perform output of the documentation into given stream.
	 */
	ostream& print(ostream& stream);

protected:
	// Need to implement the resolution function. Just call that in the base class.
	void print(ostream& stream, const TypeBase *type, unsigned int depth) {
		OutputBase::print(stream, type, depth);
	}


    void print_impl(ostream& stream, const Record *type, unsigned int depth);
    void print_impl(ostream& stream, const Array *type, unsigned int depth);
    void print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const AdHocAbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const Selection *type, unsigned int depth);
	void print_impl(ostream& stream, const Integer *type, unsigned int depth);
	void print_impl(ostream& stream, const Double *type, unsigned int depth);
	void print_impl(ostream& stream, const Bool *type, unsigned int depth);
	void print_impl(ostream& stream, const String *type, unsigned int depth);
    void print_impl(ostream& stream, const FileName *type, unsigned int depth);

private:
    /**
     * Prints value according to DefaultType
     * Respects obligatory, optional and read time flag
     *
     * @param stream Output stream
     * @param depth Depth of output
     * @param empty_val Default empty value (zero for numeric types, empty string ...)
     * @param invalid_val Flag if value is invalid for its type
     * @param has_quote Flag if value is enclosed in quotes
     */
    void print_default_value(ostream& stream, unsigned int depth, string empty_val, bool invalid_val, bool has_quote = false);

    /// temporary value of actually record type
    string key_name_;
    /// temporary value of actually reference
    string reference_;
    /// temporary value of actually record value
    Default value_;
};




/**
 * @brief Class for create documentation in Latex format using particular macros.
 *
 * @ingroup input_types
 */
class OutputLatex : public OutputBase {
public:
    OutputLatex(TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

    ostream & print(ostream& stream);

protected:

    // Need to implement the resolution function. Just call that in the base class.
    void print(ostream& stream, const TypeBase *type, unsigned int depth) {
        OutputBase::print( stream, type, depth);
    }

    void print_impl(ostream& stream, const Record *type, unsigned int depth);
    void print_impl(ostream& stream, const Array *type, unsigned int depth);
    void print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const AdHocAbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const Selection *type, unsigned int depth);
    void print_impl(ostream& stream, const Integer *type, unsigned int depth);
    void print_impl(ostream& stream, const Double *type, unsigned int depth);
    void print_impl(ostream& stream, const Bool *type, unsigned int depth);
    void print_impl(ostream& stream, const String *type, unsigned int depth);
    void print_impl(ostream& stream, const FileName *type, unsigned int depth);


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
	OutputJSONMachine(TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth)
    {
	    format_head="[\n";
	    format_tail="{}]\n";
    }

    /**
     * @brief Simple internal class for storing JSON description rewrite rule
     */
    class RewriteRule {
        public:
            /** boost regex regular expression */
            boost::regex search;

            /** simple string replacement */
            std::string replacement;

            /**
             * Constructor
             * @param _search boost::regex regular expression match pattern
             * @param _replacement std::string string _replacement
             */
            RewriteRule (boost::regex _search, std::string _replacement):
                search(_search),
                replacement(_replacement) {

            };
    };


protected:
	std::string format_hash(TypeBase::TypeHash hash);
	std::string escape_description(std::string desc);

    void print_impl(ostream& stream, const Record *type, unsigned int depth);
    void print_impl(ostream& stream, const Array *type, unsigned int depth);
    void print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const AdHocAbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const Selection *type, unsigned int depth);
    void print_impl(ostream& stream, const Integer *type, unsigned int depth);
    void print_impl(ostream& stream, const Double *type, unsigned int depth);
    void print_impl(ostream& stream, const Bool *type, unsigned int depth);
    void print_impl(ostream& stream, const String *type, unsigned int depth);
    void print_impl(ostream& stream, const FileName *type, unsigned int depth);


    /// Print all keys of AbstractRecord type or AdHocAbstractRecord type
    void print_abstract_record_keys(ostream& stream, const AbstractRecord *type, unsigned int depth);
};


/**
 * Overrides output operator for simple output of the input type tree.
 */
std::ostream& operator<<(std::ostream& stream, OutputText type_output);
std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output);
std::ostream& operator<<(std::ostream& stream, OutputLatex type_output);
std::ostream& operator<<(std::ostream& stream, OutputJSONMachine type_output);



} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_OUTPUT_HH_ */
