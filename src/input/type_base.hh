/*
 * type_base.hh
 *
 *  Created on: May 1, 2012
 *      Author: jb
 */

#ifndef TYPE_BASE_HH_
#define TYPE_BASE_HH_

#include <limits>
#include <ios>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <iomanip>

#include <boost/type_traits.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>

#include "system/global_defs.h"
#include "system/system.hh"
#include "system/exceptions.hh"
#include "system/file_path.hh"




namespace Input {

namespace Type {



using namespace std;

/**
 * Declaration of common exceptions and error info types.
 */
TYPEDEF_ERR_INFO( EI_KeyName, const string );

TYPEDEF_ERR_INFO( EI_DefaultStr, const string);
TYPEDEF_ERR_INFO( EI_TypeName, const string);
DECLARE_EXCEPTION( ExcWrongDefault, << "Default value " << EI_DefaultStr::qval
        << " do not match type: " << EI_TypeName::qval << ";\n"
        << "During declaration of the key: " << EI_KeyName::qval );




/**
 * @brief Base of classes for declaring structure of the input data.
 *
 *  Provides methods common to all types. Namely, the type name, finished status (nontrivial only for types with complex initialization - Record, AbstractRecosd, Selection)
 *  and output of the documentation.
 *
 *  @ingroup input_types
 */
class TypeBase {
public:
	typedef std::size_t TypeHash;
	typedef std::string json_string;
	typedef std::map<std::string, json_string> attribute_map;
	typedef std::pair< std::string, boost::shared_ptr<TypeBase> > ParameterPair;
	typedef std::map< std::string, TypeHash > ParameterMap;
	typedef std::pair< boost::shared_ptr<TypeBase>, ParameterMap > MakeInstanceReturnType;


    /**
     * Returns true if the type is fully specified and ready for read access. For Record and Array types
     * this say nothing about child types referenced in particular type object.
     * In particular for Record and Selection, it returns true after @p finish() method is called.
     *
     */
    virtual bool is_finished() const
    {return true;}

    /**
     * Returns true if the type is closed.
     *
     */
    virtual bool is_closed() const
    {return true;}

    /// Returns an identification of the type. Useful for error messages.
    virtual string type_name() const  { return "TypeBase"; }

    /// Returns an identification of the type in format "type_name():parent->type_name()".
    virtual string full_type_name() const  { return "TypeBase"; }

    /**
     * Returns string with Type extensive documentation. We need this to pass Type description at
     * throw points since the Type object can be deallocated during stack unrolling so it is not good idea to pass
     * pointer. Maybe we can pass smart pointers. Actually this method is used in various exceptions in json_to_storage.
     *
     * Some old note on this topic:
     *    !!! how to pass instance of descendant of TypeBase through EI -
     *  - can not pass it directly since TypeBase is not copyconstructable
     *  - can not use shared_ptr for same reason
     *  - can not use C pointers since the refered object can be temporary
     *  solutions:
     *   - consistently move TypeBase to Pimpl design
     *   - provide virtual function make_copy, that returns valid shared_ptr
     *
     */
    string desc() const;



    /**
     * Comparison of types. It compares kind of type (Integer, Double, String, Record, ..), for complex types
     * it also compares names. For arrays compare subtypes.
     */
    virtual bool operator==(const TypeBase &other) const
        { return typeid(*this) == typeid(other); }

    /// Comparison of types.
    bool operator!=(const TypeBase & other) const
        { return ! (*this == other); }

    /**
     *  Destructor
     */
    virtual ~TypeBase();



    /// Finishes all types registered in type repositories.
    static void lazy_finish();


    /**
     * Finish method. Finalize construction of "Lazy types": Record, Selection, AbstractRecord, and Array.
     * These input types are typically defined by means
     * of static variables, whose order of initialization is not known a priori. Since e.g. a Record can link to other
     * input types through its keys, these input types cannot be accessed directly at the initialization phase.
     * The remaining part of initialization can be done later, typically from main(), by calling the method finish().
     *
     * Finish try to convert all raw pointers pointing to lazy types into smart pointers to valid objects. If there
     * are still raw pointers to not constructed objects the method returns false.
     */
    virtual bool finish()
    { return true; };

    /**
     * For types that can be initialized from a default string, this method check
     * validity of the default string. For invalid string an exception is thrown.
     *
     * Return false if the validity can not be decided due to presence of unconstructed types (Record, Selection)
     */
    virtual bool valid_default(const string &str) const =0;

    /**
     * Hash of the type specification. Provides unique id computed from its
     * content (definition) so that same types have same hash.
     */
    virtual TypeHash content_hash() const =0;

    /// Add attribute to map
    void add_attribute(std::string name, json_string val);

    /// Print JSON output of attributes to @p stream.
    void write_attributes(ostream& stream) const;

    /// Create instance of generic type
    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const =0;

    /// Return if type is generic.
    inline bool is_generic() const
    { return generic_; };

protected:

    /**
     * Default constructor.
     */
    TypeBase();

    /**
     * Copy constructor.
     */
    TypeBase(const TypeBase& other);




    /**
     * Type of hash values used in associative array that translates key names to indices in Record and Selection.
     *
     * For simplicity, we currently use whole strings as "hash".
     */
    typedef string KeyHash;

    /// Hash function.
    inline static KeyHash key_hash(const string &str) {
        return (str);
    }

    /**
     * Check that a @p key is valid identifier, i.e. consists only of valid characters, that are lower-case letters, digits and underscore,
     * we allow identifiers starting with a digit, but it is discouraged since it slows down parsing of the input file.
     */
    static bool is_valid_identifier(const string& key);

    /// Check if JSON string is valid
    bool validate_json(json_string str) const;

    /**
     * Add attributes to hash of the type specification.
     *
     * Method must be called in content_hash() method of TypeBase descendants.
     */
    void attribute_content_hash(std::size_t &seed) const;

    /// Merge values @p other map to @p parameter_map_
    void add_to_parameter_map(ParameterMap other);

    /// Create JSON output from @p parameter_map_ formatted as attribute.
    json_string print_parameter_map_to_json();

    /// map of type attributes (e. g. input_type, name etc.)
    boost::shared_ptr<attribute_map> attributes_;

    /**
     * Map of used parameters in generic type (for non-generic types is empty).
     *
     * Contains all parameter pairs used in type and in all descendants (in subtree of type).
     */
    mutable ParameterMap parameter_map_;

    /// Flag determining if type is generic (e. g. Record is generic if has key of type Parameter).
    bool generic_;

    friend class Array;
    friend class Record;
};

/**
 * For convenience we provide also redirection operator for output documentation of Input:Type classes.
 */
std::ostream& operator<<(std::ostream& stream, const TypeBase& type);


class Record;
class Selection;


/**
 * @brief Class for declaration of inputs sequences.
 *
 * The type is fully specified after its constructor is called. All elements of the Array has same type, however you
 * can use elements of AbstractRecord.
 *
 * If you not disallow Array size 1, the input reader will try to convert any other type
 * on input into array with one element, e.g.
 @code
     int_array=1        # is equivalent to
     int_array=[ 1 ]
 @endcode
 *
 * @ingroup input_types
 */
class Array : public TypeBase {
	friend class OutputBase;

protected:

    class ArrayData  {
    public:

    	ArrayData(unsigned int min_size, unsigned int max_size)
    	: lower_bound_(min_size), upper_bound_(max_size), finished(false)
    	{}

    	ArrayData(const ArrayData &other)
    	: type_of_values_(other.type_of_values_), lower_bound_(other.lower_bound_),
		  upper_bound_(other.upper_bound_), finished(false)
    	{}

    	bool finish();

    	boost::shared_ptr<TypeBase> type_of_values_;
    	unsigned int lower_bound_, upper_bound_;
    	bool finished;

    };

public:
    /**
     * Constructor with a @p type of array items given as pure reference. In this case \p type has to by descendant of \p TypeBase different from
     * 'complex' types @p Record and @p Selection. You can also specify minimum and maximum size of the array.
     *
     * Possibly modification: If we need set attributes to Array, we can't need TypeBase::add_attribute method (requires non-closed type).
     * Solution: We add to constructor optional parameter of type vector< pair<std::string, json_string> > what allows set all attributes.
     */
    template <class ValueType>
    Array(const ValueType &type, unsigned int min_size=0, unsigned int max_size=std::numeric_limits<unsigned int>::max() );

    TypeHash content_hash() const override;

    /// Finishes initialization of the Array type because of lazy evaluation of type_of_values.
    virtual bool finish();

    virtual bool is_finished() const {
        return data_->finished; }

    /// Getter for the type of array items.
    inline const TypeBase &get_sub_type() const {
        return *data_->type_of_values_; }

    /// Checks size of particular array.
    inline bool match_size(unsigned int size) const {
        return size >=data_->lower_bound_ && size<=data_->upper_bound_; }

    /// @brief Implements @p Type::TypeBase::type_name. Name has form \p array_of_'subtype name'
    virtual string type_name() const;

    /// @brief Implements @p Type::TypeBase::full_type_name.
    virtual string full_type_name() const;

    /// @brief Implements @p Type::TypeBase::operator== Compares also subtypes.
    virtual bool operator==(const TypeBase &other) const;

    /**
     *  Default values for an array creates array containing one element
     *  that is initialized by given default value. So this method check
     *  if the default value is valid for the sub type of the array.
     */
    virtual bool valid_default(const string &str) const;

    // Implements @p TypeBase::make_instance.
    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const;

    /// Create deep copy of Array (copy all data stored in shared pointers etc.)
    Array deep_copy() const;

protected:

    /// Handle to the actual array data.
    boost::shared_ptr<ArrayData> data_;
private:
    /// Forbids default constructor in order to prevent empty data_.
    Array();
};



/**
 * @brief Base of all scalar types.
 *
 *  @ingroup input_types
 */
class Scalar : public TypeBase {
public:

	virtual string full_type_name() const;

};


/**
 * @brief Class for declaration of the input of type Bool.
 *
 * String names of boolean values are \p 'true' and \p 'false'.
 *
 * @ingroup input_types
 */
class Bool : public Scalar {
public:
    Bool()
	{}

    TypeHash content_hash() const   override;


    bool from_default(const string &str) const;

    virtual string type_name() const;

    virtual bool valid_default(const string &str) const;

    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const;
};


/**
 * @brief Class for declaration of the integral input data.
 *
 * The data are stored in an \p signed \p int variable. You can specify bounds for the valid input data.
 *
 * @ingroup input_types
 */
class Integer : public Scalar {
	friend class OutputBase;

public:
    Integer(int lower_bound=std::numeric_limits<int>::min(), int upper_bound=std::numeric_limits<int>::max())
	: lower_bound_(lower_bound), upper_bound_(upper_bound)
	{}

    TypeHash content_hash() const   override;

    /**
     * Returns true if the given integer value conforms to the Type::Integer bounds.
     */
    bool match(std::int64_t value) const;

    /**
     * As before but also returns converted integer in @p value.
     */
    int from_default(const string &str) const;
    /// Implements  @p Type::TypeBase::valid_defaults.
    virtual bool valid_default(const string &str) const;

    virtual string type_name() const;

    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const;
private:

    std::int64_t lower_bound_, upper_bound_;

};


/**
 * @brief Class for declaration of the input data that are floating point numbers.
 *
 * The data are stored in an \p double variable. You can specify bounds for the valid input data.
 *
 * @ingroup input_types
 */
class Double : public Scalar {
	friend class OutputBase;

public:
    Double(double lower_bound= -std::numeric_limits<double>::max(), double upper_bound=std::numeric_limits<double>::max())
	: lower_bound_(lower_bound), upper_bound_(upper_bound)
	{}

    TypeHash content_hash() const   override;

    /**
     * Returns true if the given integer value conforms to the Type::Double bounds.
     */
    bool match(double value) const;

    /// Implements  @p Type::TypeBase::valid_defaults.
    virtual bool valid_default(const string &str) const;

    /**
     * As before but also returns converted integer in @p value.
     */
    double from_default(const string &str) const;

    virtual string type_name() const;

    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const;
private:


    double lower_bound_, upper_bound_;

};



/**
 * Just for consistency, but is essentially same as Scalar.
 *
 * @ingroup input_types
 */
class String : public Scalar {
public:
    virtual string type_name() const;

    TypeHash content_hash() const   override;


    string from_default(const string &str) const;

    /**
     * Particular descendants can check validity of the string.
     */
    virtual bool match(const string &value) const;

    /// Implements  @p Type::TypeBase::valid_defaults.
    virtual bool valid_default(const string &str) const;

    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const;
};


/**
 * @brief Class for declaration of the input data that are file names.
 *
 * We strictly distinguish filenames for input and output files.
 *
 * @ingroup input_types
 */
class FileName : public String {
public:

	TypeHash content_hash() const   override;

    /**
     * Factory function for declaring type FileName for input files.
     */
    static FileName input()
    { return FileName(::FilePath::input_file); }

    /**
     * Factory function for declaring type FileName for input files.
     */
    static FileName output()
    { return FileName(::FilePath::output_file); }

    virtual string type_name() const;

    virtual bool operator==(const TypeBase &other) const
    { return  typeid(*this) == typeid(other) &&
                     (type_== static_cast<const FileName *>(&other)->get_file_type() );
    }

    /// Checks relative output paths.
    virtual bool match(const string &str) const;


    /**
     * Returns type of the file input/output.
     */
    ::FilePath::FileType get_file_type() const {
        return type_;
    }



private:
    ::FilePath::FileType    type_;

    /// Forbids default constructor.
    FileName() {}

    /// Forbids direct construction.
    FileName(enum ::FilePath::FileType type)
    : type_(type)
    {}

};

} // closing namespace Type
} // closing namespace Input





#endif /* TYPE_BASE_HH_ */
