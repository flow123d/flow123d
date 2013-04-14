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
/***
 * TODO:
 *
 * oprava get_xy_data

simplify implementation of has_type a souvisejicich metod, odstranit inline
metody pro pristup k ProcessedData presunout z OutputBase do ProcessedData
metody:
 bool was_written( string full_name)
 void mark_written( string full_name)
... pripadne spojit
aplikuje filter_ na full_name :
    boost::regex_replace(str_out.begin(), str_in.begin(), str_in.end(), filter_, "$&", ...)

v OutputBase metoda OutputBase & set_filter(const string &reg_exp) -> nastavuje reg_exp v ProcessedTypes

dalsi mapa full_type_name - moznost filtorvat vystupy pomoci regularniho vyrazu

/// implement both DFS and BFS print, add methods push() and pop() to OutputBase

 *
 */

/**
 * Base abstract class for output description of the Input::Type tree.
 * Output into various formats is implemented by derived classes.
 *
 * Usage:
 * cout << OutputText( &my_record, 3) << "konec" << endl;
 *
 */
class OutputBase {
public:

    /**
     * Performs output of the documentation into given @p stream.
     * Returns reference to the same stream.
     */
    virtual ostream& print(ostream& stream);

	/// Initialize filter_; alokace filter_
	void set_filter(string regex_filter);

protected:
    /// Types of documentation output
    enum DocumentationType {
        key_record,
        full_record
    };

    /**
     * Structure for flags about output of one TypeBase object in input tree
     */
    struct Key {
    	unsigned int key_index;          	///< Position inside the record.
    	const void * type_;              	///< Pointer to type.
    	mutable bool extensive_doc_;     	///< Flag captures if extensive documentation of type was printed.
    	mutable string reference_;       	///< Reference to type.
    };

    /**
     * Public typedef of constant iterator into array of keys.
     */
    typedef std::vector<struct Key>::const_iterator KeyIter;


    /// Padding of new level of printout, used where we use indentation.
    static const unsigned int padding_size = 4;


    /**
     * Constructor
     *
     * @param type Stores input sequence
     * @param depth Depth of output
     */
    OutputBase(const TypeBase *type, unsigned int depth = 0);


    // destructor
    virtual ~OutputBase();

    // data getters
    void get_array_sizes(Array array, unsigned int &lower , unsigned int &upper );
    void get_record_key(Record rec, unsigned int key_idx, Record::Key &key);
    void get_integer_bounds(Integer integer, int &lower , int &upper );
    void get_double_bounds(Double dbl, double &lower , double &upper );
    void get_parent_ptr(Record rec, boost::shared_ptr<AbstractRecord> &parent_ptr);
    void get_array_type(Array array, boost::shared_ptr<const TypeBase> &arr_type);
    const void * get_record_data(const Record *rec);
    const void * get_abstract_record_data(const AbstractRecord *a_rec);
    const void * get_selection_data(const Selection *sel);
    const void * get_array_data(const Array *array);
    const void * get_type_base_data(const TypeBase *type);


    /**
     * Perform resolution according to actual @p type (using typeid) and call particular print_impl method.
     */
    void print(ostream& stream, const TypeBase *type, unsigned int depth);


    /**
     *  following methods realize output in particular format
     *  using getters from the base class OutputBase
     */
    virtual void print_impl(ostream& stream, const Record *type, unsigned int depth) = 0;
    virtual void print_impl(ostream& stream, const Array *type, unsigned int depth) = 0;
    virtual void print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth) = 0;
    virtual void print_impl(ostream& stream, const Selection *type, unsigned int depth) = 0;
	virtual void print_impl(ostream& stream, const Integer *type, unsigned int depth) = 0;
	virtual void print_impl(ostream& stream, const Double *type, unsigned int depth) = 0;
	virtual void print_impl(ostream& stream, const Bool *type, unsigned int depth) = 0;
	virtual void print_impl(ostream& stream, const String *type, unsigned int depth) = 0;
    virtual void print_impl(ostream& stream, const FileName *type, unsigned int depth) = 0;

    /**
     * Write out a string with given padding of every new line.
     *
     * @param stream Output stream
     * @param str Printed description
     * @param hash_count Count of '#' chars in description
     */
    //virtual void write_description(std::ostream& stream, const string& str, unsigned int hash_count = 1) = 0;
    /**
     * Output indented multi-line string.
     */
    void write_description(std::ostream& stream, const string& str, unsigned int padding, unsigned int hash_count = 1);


    /**
     * Returns true if the ProcessedTypes contains key of given type and key has true flag extensive_doc_.
     */
    bool has_type_extensive(const void * type) const;

    /**
     * Returns reference_ string of key of given type.
     */
    const string get_reference(const void * type) const;

    /**
     * Write value stored in dft.
     *
     * Enclose value in quotes if it's needed or write info that value is optional or obligatory.
     */
    void write_value(std::ostream& stream, Default dft);



    /// Object for which is created printout
    const TypeBase *type_;
    /// Depth of printout
    unsigned int depth_;
    /// Type of documentation output
    DocumentationType doc_type_;

    /// temporary value for printout of description (used in std::setw function)
    unsigned int size_setw_;

    /**
     * Internal data class
     */
    class ProcessedTypes {
    public:

    	/**
    	 * Declare a processed type and its flags.
    	 *
    	 * Pointer to type must be unique in map. If pointer exists type is not added and method returns false.
    	 */
    	bool add_type(const void *type, bool extensive_doc, string reference);

    	/// Declare a processed type with default values of flags
    	bool add_type(const void *type);

    	/// Clear all data of processed types
    	void clear();

        /**
         * Interface to mapping key -> index. Returns index (in continuous array) for given type.
         */
        unsigned int type_index(const void * type) const;

        /**
         * Remove key of given type.
         */
        void remove_type(const void * type);

        /**
         * Set reference_ string of key of given type.
         */
        void set_reference(const void * type, const string& ref);

        /**
         * Set value to extensive_doc_ flag of key of given type.
         */
        void set_extensive_flag(const void * type, bool val = true);

    	/// Deallocate filter_ if it was allocated.
    	~ProcessedTypes();

        /**
         * Returns true if the ProcessedTypes contains type of given full_name
         * in full_type_names set and regular expression filter is initialized.
         */
    	bool was_written(string full_name);

    	/**
    	 * Marks full_name of type as written.
    	 * Regular expression filter must be initialized!
    	 */
    	void mark_written(string full_name);

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
    };

    /// Stores flags and references of processed type
    ProcessedTypes doc_flags_;

};


/**********************************************************************************************************************/

/**
 * Class for create text documentation
 */
class OutputText : public OutputBase {
public:
	OutputText(const TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

protected:

    void print_impl(ostream& stream, const Record *type, unsigned int depth);
    void print_impl(ostream& stream, const Array *type, unsigned int depth);
    void print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const Selection *type, unsigned int depth);
	void print_impl(ostream& stream, const Integer *type, unsigned int depth);
	void print_impl(ostream& stream, const Double *type, unsigned int depth);
	void print_impl(ostream& stream, const Bool *type, unsigned int depth);
	void print_impl(ostream& stream, const String *type, unsigned int depth);
    void print_impl(ostream& stream, const FileName *type, unsigned int depth);


};







/**
 * Class for create and JSON template documentation
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
    void print_impl(ostream& stream, const Selection *type, unsigned int depth);
	void print_impl(ostream& stream, const Integer *type, unsigned int depth);
	void print_impl(ostream& stream, const Double *type, unsigned int depth);
	void print_impl(ostream& stream, const Bool *type, unsigned int depth);
	void print_impl(ostream& stream, const String *type, unsigned int depth);
    void print_impl(ostream& stream, const FileName *type, unsigned int depth);

    //void write_description(std::ostream& stream, const string& str, unsigned int hash_count = 1);

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
    /// temporary value of actually record value
    Default value_;
};




/**
 * Class for create and Latex documentation
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
    void print_impl(ostream& stream, const Selection *type, unsigned int depth);
    void print_impl(ostream& stream, const Integer *type, unsigned int depth);
    void print_impl(ostream& stream, const Double *type, unsigned int depth);
    void print_impl(ostream& stream, const Bool *type, unsigned int depth);
    void print_impl(ostream& stream, const String *type, unsigned int depth);
    void print_impl(ostream& stream, const FileName *type, unsigned int depth);


};


/**
 * Overrides output operator for simple output of the input type tree.
 */
std::ostream& operator<<(std::ostream& stream, OutputText type_output);
std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output);
std::ostream& operator<<(std::ostream& stream, OutputLatex type_output);



} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_OUTPUT_HH_ */
