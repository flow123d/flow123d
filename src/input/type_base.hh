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
 * @file    type_base.hh
 * @brief   
 */

#ifndef TYPE_BASE_HH_
#define TYPE_BASE_HH_

#include <limits>
#include <ios>
//#include <set>
#include <map>
#include <vector>
#include <string>
//#include <iomanip>

#include <stdint.h>                  // for int64_t
#include <boost/exception/info.hpp>  // for error_info::~error_info<Tag, T>
#include <memory>                    // for shared_ptr, __shared_ptr
#include <utility>                   // for pair


//#include <boost/type_traits.hpp>
//#include <boost/tokenizer.hpp>
//#include <boost/algorithm/string.hpp>

//#include "system/global_defs.h"
//#include "system/system.hh"
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
TYPEDEF_ERR_INFO( EI_Desc, const string);
DECLARE_EXCEPTION( ExcWrongDefaultJSON, << "Consistency Error: Not valid JSON of Default value "
		<< EI_DefaultStr::qval << " of type " << EI_TypeName::qval << ";\n"
		<< "During declaration of the key: " << EI_KeyName::qval );
DECLARE_EXCEPTION( ExcWrongDefault, << "Consistency Error: " << EI_Desc::val << "Default value "
		<< EI_DefaultStr::qval << " do not match type: " << EI_TypeName::qval << ";\n"
        << "During declaration of the key: " << EI_KeyName::qval );
DECLARE_EXCEPTION( ExcUnknownDescendant, << "Unknown descendant of TypeBase class, name: " << EI_TypeName::qval );




/**
 * FinishStatus manages finish of elements in IST.
 *
 * Specifies:
 *  - finish status of IST elements
 *  - type of executing finish
 */
enum FinishStatus {
	none_,        //< unfinished element / this type can't be used as type of finish
	in_perform_,  //< finish of element performs (used in recursion check) / this type can't be used as type of finish
	regular_,     //< finished element of IST / finish of IST executed recursively from root element
	generic_,     //< finished element in generic subtree / finish of generic subtree (executed from Instance object)
	delete_       //< finished element marked as deleted (that doesn't appear in IST) / finish of unused elements (that can be removed)
};


/**
 * @brief Base of classes for declaring structure of the input data.
 *
 *  Provides methods and class members common to all types. Namely, the type name, finished status (nontrivial only
 *  for types with complex initialization - Record, AbstractRecosd, Selection), attributes, data of generic types
 *  and output of the documentation.
 *
 *  @ingroup input_types
 */
class TypeBase {
public:
	/// Type returned by content_hash methods.
	typedef std::size_t TypeHash;
	const static TypeHash none_hash = std::numeric_limits<TypeHash>::max();

	/// String stored in JSON format.
	typedef std::string json_string;
	/// Defines map of Input::Type attributes.
	typedef std::map<std::string, json_string> attribute_map;

	/// Defines pairs of (name, Input::Type), that are used for replace of parameters in generic types.
	typedef std::pair< std::string, std::shared_ptr<TypeBase> > ParameterPair;
	/// Define vector of parameters passed to the overloaded make_instance method.
	typedef std::vector<ParameterPair> ParameterVector;

	/// Defines map of used parameters
	typedef std::map< std::string, TypeHash > ParameterMap;
	/// Return type of make_instance methods, contains instance of generic type and map of used parameters
	typedef std::pair< std::shared_ptr<TypeBase>, ParameterMap > MakeInstanceReturnType;


    /**
     * @brief Returns true if the type is fully specified and ready for read access.
     *
     * For Record and Array types this say nothing about child types referenced in particular type object.
     * In particular for Record and Selection, it returns true after @p finish() method is called.
     */
    virtual FinishStatus finish_status() const;

    virtual bool is_finished() const;

    /**
     * @brief Returns true if the type is closed.
     */
    virtual bool is_closed() const;

    /// Returns an identification of the type. Useful for error messages.
    virtual string type_name() const;
    /// Returns an identification of the class. Useful for output of the documentation.
    virtual string class_name() const;

    /**
     * @brief Returns string with Type extensive documentation.
     *
     * We need this to pass Type description at throw points since the Type object can be deallocated during stack
     * unrolling so it is not good idea to pass pointer. Maybe we can pass smart pointers. Actually this method is
     * used in various exceptions in json_to_storage.
     *
     * TODO: Some old note on this topic:
     *    !!! how to pass instance of descendant of TypeBase through EI -
     *  - can not pass it directly since TypeBase is not copyconstructable
     *  - can not use shared_ptr for same reason
     *  - can not use C pointers since the referred object can be temporary
     *  solutions:
     *   - consistently move TypeBase to Pimpl design
     *   - provide virtual function make_copy, that returns valid shared_ptr
     *
     */
    string desc() const;



    /**
     * @brief Comparison of types.
     *
     * It compares kind of type (Integer, Double, String, Record, ..), for complex types it also compares names.
     * For arrays compare subtypes.
     */
    virtual bool operator==(const TypeBase &other) const;

    /// Comparison of types.
    bool operator!=(const TypeBase & other) const;

    ///  Destructor
    virtual ~TypeBase();



    /**
     * @brief Finishes and marks all types registered in type repositories and unused in IST.
     *
     * Steps of this method:
     *  1) finishes all types unused in IST, marks them as deleted (FinishStatus::deleted_)
     *  2) iterates these deleted types once more, checks shared pointers link to Record keys, descendants of Abstract etc.
     *     - if count == 1, it's OK, method resets shared pointer
     *     - in other cases throws error
     */
    static void delete_unfinished_types();


    /**
     * @brief Finish method. Finalize construction of "Lazy types": Record, Selection, Abstract and generic type.
     *
     * These input types are typically defined by means of static generating methods, whose allows initialization
     * any of other input types. Since e.g. a Record can link to other input types through its keys at the
     * initialization phase. But some setting (e.g. add descendants to Abstract) can't be done in initialization
     * phase. Therefore, the remaining part of initialization can be done later, in finalization phase, typically
     * from main(), by calling the method finish().
     *
     * Method allows finish in different cases that is managed by @p finish_type:
     *
     * - base finish of IST requires FinishStatus::regular_ and typically finish of root type is called in this case
     * - finish of generic types can be different of other Input::Types (e. g. for Record) and needs set @p finish_type
     *   to FinishStatus::generic_
     * - input types unused in IST can be finished with parameter FinishStatus::delete_. These types can be then safely
     *   deleted from \p Input::TypeRepository (see also \p delete_unfinished_types)
     *
     */
    virtual FinishStatus finish(FinishStatus finish_type = FinishStatus::regular_);

    /**
     * @brief Hash of the type specification.
     *
     * Provides unique id computed from its content (definition) so that same types have same hash.
     *
     * Hash is counted using type, name and other class members specific for descendants.
     */
    virtual TypeHash content_hash() const =0;

    /**
     * @brief Format given hash for output.
     *
     * Use hex format and double quotas.
     */
    static std::string hash_str(TypeHash hash);

    /// Format the hash of this type.
    inline std::string hash_str() const {
        return hash_str(content_hash());
    }

    /**
     * Create instance of generic type.
     *
     * Replace parameters in input tree by type stored in @p vec.
     */
    virtual MakeInstanceReturnType make_instance(ParameterVector vec = ParameterVector()) =0;

    /// Indicates if type is marked with flag @p root_of_generic_subtree_
    inline bool is_root_of_generic_subtree() const {
    	return root_of_generic_subtree_;
    }

protected:

    /// The default constructor.
    TypeBase();

    /// Copy constructor.
    TypeBase(const TypeBase& other);

    /**
     * @brief Add attribute of given @p name to attribute map.
     *
     * Parameter @p val must contain valid JSON string.
     */
    void add_attribute_(std::string name, json_string val);



    /**
     * @brief The type of hash values used in associative array that translates key names to indices in Record and Selection.
     *
     * For simplicity, we currently use whole strings as "hash".
     */
    typedef string KeyHash;

    /// Hash function.
    inline static KeyHash key_hash(const string &str) {
        return (str);
    }

    /**
     * @brief Check that a @p key is valid identifier.
     *
     * I.e. consists only of valid characters, that are lower-case letters, digits and underscore,
     * we allow identifiers starting with a digit, but it is discouraged since it slows down parsing of the input file.
     */
    static bool is_valid_identifier(const string& key);

    /// Check if @p str is valid JSON string
    bool validate_json(json_string str) const;

    /// Create JSON output from @p parameter_map formatted as value of attribute.
    json_string print_parameter_map_to_json(ParameterMap parameter_map) const;

    /**
     * JSON format of the vector of parameter names. Used for generic composed types.
     * TODO: use ParameterVector instead of ParameterMap and merge this method with the previous one.
     */
    json_string print_parameter_map_keys_to_json(ParameterMap param_map) const;

    /**
     * Extract and set attributes from the ParameterMap containing parameters of the subtree.
     * Used in implementations of the make_instance method.
     */
    void set_generic_attributes(ParameterMap param_map);

    /**
     * Copy attributes to the instance of the generic type. Remove all internal attributes starting with '_'.
     */
    void copy_attributes(attribute_map other_attributes);


    /// map of type attributes (e. g. input_type, name etc.)
    std::shared_ptr<attribute_map> attributes_;

    /// flag is true if type should be root of generic subtree
    bool root_of_generic_subtree_;

    /// hash string of generic type if type is derived, or empty string
    TypeHash generic_type_hash_;

    /// map of parameters if type is part of generic subtree
    ParameterMap parameter_map_;

    friend class Array;
    friend class Record;
    friend class OutputBase;
};

/**
 * @brief For convenience we provide also redirection operator for output documentation of Input:Type classes.
 */
std::ostream& operator<<(std::ostream& stream, const TypeBase& type);


//class Record;
//class Selection;


/**
 * @brief Class for declaration of inputs sequences.
 *
 * The type is fully specified after its constructor is called. All elements of the Array has same type, however you
 * can use elements of Abstract.
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

	/**
	 * @brief Actual data of the Array.
	 */
    class ArrayData  {
    public:

    	/// Constructor
    	ArrayData(unsigned int min_size, unsigned int max_size);
    	/// Finishes initialization of the ArrayData.
    	FinishStatus finish(FinishStatus finish_type = FinishStatus::regular_);
    	/// Type of Array
    	std::shared_ptr<TypeBase> type_of_values_;
    	/// Minimal size of Array
    	unsigned int lower_bound_;
    	/// Maximal size of Array
    	unsigned int upper_bound_;
    	/// Flag specified if Array is finished
    	FinishStatus finish_status;

    };

public:
    /**
     * @brief Constructor with a @p type of array items given as pure reference.
     *
     * In this case \p type has to by descendant of \p TypeBase different from 'complex' types
     * @p Record and @p Selection. You can also specify minimum and maximum size of the array.
     */
    template <class ValueType>
    Array(const ValueType &type, unsigned int min_size=0, unsigned int max_size=std::numeric_limits<unsigned int>::max() );

    /**
     * @brief Constructor with a shared pointer @p type of array.
     */
    Array(std::shared_ptr<TypeBase> type, unsigned int min_size=0, unsigned int max_size=std::numeric_limits<unsigned int>::max() );

    /**
     * @brief Implements @p TypeBase::content_hash.
     *
     * Hash is calculated by type name, bounds, hash of stored type and hash of attributes.
     */
    TypeHash content_hash() const override;

    /// Finishes initialization of the Array type because of lazy evaluation of type_of_values.
    FinishStatus finish(FinishStatus finish_type = FinishStatus::regular_) override;

    /// Override @p Type::TypeBase::finish_status.
    FinishStatus finish_status() const override;

    /// Override @p Type::TypeBase::is_finished.
    bool is_finished() const override;

    /// Getter for the type of array items.
    inline const TypeBase &get_sub_type() const {
        return *data_->type_of_values_; }

    /// Checks size of particular array.
    inline bool match_size(unsigned int size) const {
        return size >=data_->lower_bound_ && size<=data_->upper_bound_; }

    /**
     * @brief Implements @p Type::TypeBase::type_name.
     *
     * Name has form \p array_of_'subtype_name'.
     */
    string type_name() const override;
    /// Override @p Type::TypeBase::class_name.
    string class_name() const override;

    /// @brief Implements @p Type::TypeBase::operator== Compares also subtypes.
    bool operator==(const TypeBase &other) const override;

    /// Implements @p TypeBase::make_instance.
    MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) override;

    /**
     * @brief Create deep copy of Array.
     *
     * Copy all data stored in shared pointers etc.
     */
    Array deep_copy() const;

protected:

    /// Handle to the actual array data.
    std::shared_ptr<ArrayData> data_;
private:
    /// Forbids default constructor in order to prevent empty data_.
    Array();
};



/**
 * @brief Base of all scalar types.
 *
 *  @ingroup input_types
 */
class Scalar : public TypeBase {};


/**
 * @brief Class for declaration of the input of type Bool.
 *
 * String names of boolean values are \p 'true' and \p 'false'.
 *
 * @ingroup input_types
 */
class Bool : public Scalar {
public:
	/// Constructor.
    Bool()
	{}

    /// Implements @p TypeBase::content_hash.
    TypeHash content_hash() const   override;


    /**
     * @brief Implements @p Type::TypeBase::type_name.
     *
     * Name has form \p Bool.
     */
    string type_name() const override;
    /// Override @p Type::TypeBase::class_name.
    string class_name() const override;

    /// Implements @p TypeBase::make_instance.
    MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) override;
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
	/**
	 * @brief Constructor.
	 *
	 * You can also specify minimum and maximum value.
	 */
    Integer(int lower_bound=std::numeric_limits<int>::min(), int upper_bound=std::numeric_limits<int>::max());

    /**
     * @brief Implements @p TypeBase::content_hash.
     *
     * Hash is calculated by type name and bounds.
     */
    TypeHash content_hash() const   override;

    /**
     * @brief Check valid value of Integer.
     *
     * Returns true if the given integer value conforms to the Type::Integer bounds.
     */
    bool match(std::int64_t value) const;

    /**
     * @brief Implements @p Type::TypeBase::type_name.
     *
     * Name has form \p Integer.
     */
    string type_name() const override;
    /// Override @p Type::TypeBase::class_name.
    string class_name() const override;

    /// Implements @p TypeBase::make_instance.
    MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>())  override;
private:

    std::int64_t lower_bound_; ///< Minimal value of Integer.
    std::int64_t upper_bound_; ///< Maximal value of Integer.

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
	/**
	 * @brief Constructor.
	 *
	 * You can also specify minimum and maximum value.
	 */
    Double(double lower_bound= -std::numeric_limits<double>::max(), double upper_bound=std::numeric_limits<double>::max());

    /**
     * @brief Implements @p TypeBase::content_hash.
     *
     * Hash is calculated by type name and bounds.
     */
    TypeHash content_hash() const   override;

    /// Returns true if the given integer value conforms to the Type::Double bounds.
    bool match(double value) const;

    /**
     * @brief Implements @p Type::TypeBase::type_name.
     *
     * Name has form \p Double.
     */
    string type_name() const override;
    /// Override @p Type::TypeBase::class_name.
    string class_name() const override;

    /// Implements @p TypeBase::make_instance.
    MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) override;
private:

    double lower_bound_; ///< Minimal value of Integer.
    double upper_bound_; ///< Maximal value of Integer.

};



/**
 * @brief Class for declaration of the input data that are in string format.
 *
 * Just for consistency, but is essentially same as Scalar.
 *
 * @ingroup input_types
 */
class String : public Scalar {
public:
    /**
     * @brief Implements @p Type::TypeBase::type_name.
     *
     * Name has form \p String.
     */
    virtual string type_name() const override;
    /// Override @p Type::TypeBase::class_name.
    string class_name() const override;

    /// Implements @p TypeBase::content_hash.
    TypeHash content_hash() const   override;

    /// Particular descendants can check validity of the string.
    /// This returns always true.
    virtual bool match(const string &value) const;

    /// Implements @p TypeBase::make_instance.
    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) override;
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

    /// Implements @p TypeBase::content_hash.
	TypeHash content_hash() const   override;

    /**
     * @brief The factory function for declaring type FileName for input files.
     */
    static FileName input();

    /**
     * @brief The factory function for declaring type FileName for input files.
     */
    static FileName output();

    /**
     * @brief Implements @p Type::TypeBase::type_name.
     *
     * Name has form \p FileName_input or \p FileName_output
     */
    string type_name() const override;
    /// Override @p Type::TypeBase::class_name.
    string class_name() const override;

    /// Comparison of types.
    bool operator==(const TypeBase &other) const override;

    /// Checks relative output paths.
    bool match(const string &str) const override;


    /**
     * @brief Returns type of the file input/output.
     */
    ::FilePath::FileType get_file_type() const;

    /// Implements @p TypeBase::make_instance.
    MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) override;



private:
    /// The type of file (input or output).
    ::FilePath::FileType    type_;

    /// Forbids default constructor.
    FileName();

    /// Forbids direct construction.
    FileName(enum ::FilePath::FileType type);

};

} // closing namespace Type
} // closing namespace Input





#endif /* TYPE_BASE_HH_ */
