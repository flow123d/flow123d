/*
 * type_output.hh
 *
 *  Created on: Nov 26, 2012
 *      Author: jb
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
    /// Gets range of integer
    void get_integer_bounds(Integer integer, int &lower , int &upper );
    /// Gets range of double
    void get_double_bounds(Double dbl, double &lower , double &upper );
    /// Gets pointer of parent AbstractRecord for given Record
    void get_parent_vec(Record rec, std::vector< boost::shared_ptr<AbstractRecord> > &parent_vec);
    /// Gets pointer of inner type for given Array
    void get_array_type(Array array, boost::shared_ptr<TypeBase> &arr_type);
    /// Gets values of default for given record key
    void get_default(Record::KeyIter it, string &type, string &value);
    /// Gets description of the given selection type.
    const string & get_selection_description(const Selection *sel);
    /// Gets parent_name_ of the given AdHocAbstractRecord type.
    const string & get_adhoc_parent_name(const AdHocAbstractRecord *a_rec);


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

    	/// Clear all data of processed types
    	void clear();


    	/// Destructor
    	~ProcessedTypes();


        /**
         * Returns true if the type was printed out
         *
         * Checks if the ProcessedTypes contains key of given type and key has true flag extensive_doc_.
         */
        bool was_written(std::size_t hash)
        {
            bool in_set = ( output_hash.find(hash) != output_hash.end() );
            if (! in_set) output_hash.insert(hash);
            return in_set;
        }


        /// Set of hashes of outputed types. Should replace keys.
        std::set<std::size_t> output_hash;
    };

    /// Stores flags and references of processed type
    ProcessedTypes doc_flags_;

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
std::ostream& operator<<(std::ostream& stream, OutputJSONMachine type_output);



} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_OUTPUT_HH_ */
