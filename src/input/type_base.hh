/*
 * type_base.hh
 *
 *  Created on: May 1, 2012
 *      Author: jb
 */

#ifndef TYPE_BASE_HH_
#define TYPE_BASE_HH_


#include "system/system.hh"
#include "system/exceptions.hh"
#include "system/file_path.hh"


#include <limits>
#include <ios>
#include <map>
#include <vector>
#include <string>
#include <iomanip>

#include <boost/type_traits.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>




namespace Input {

namespace Type {



using namespace std;



TYPEDEF_ERR_INFO( EI_KeyName, const string );

TYPEDEF_ERR_INFO( EI_DefaultStr, const string);
TYPEDEF_ERR_INFO( EI_TypeName, const string);
DECLARE_EXCEPTION( ExcWrongDefault, << "Default value " << EI_DefaultStr::qval
        << " do not match type: " << EI_TypeName::qval << ";\n"
        << "During declaration of the key: " << EI_KeyName::qval );


/**
 * @brief Default specifies default value of keys of a @p Record.
 *
 * It contains type of default value and possibly the value itself as a string.
 *
 * We prefer to use only default values specified at declaration since only those can by documented as part of
 * Record type specification.
 *
 * @ingroup input
 */
class Default {
private:
    /**
     * Possible types of default values.
     */
    enum DefaultType {
        declaration,        ///< Default value given at declaration time.
        optional_type,      ///< No default value, optional key. This is default type of the Default.
        obligatory_type     ///< No default value, obligatory key.
    };
public:


    /**
     * Constructor with given default value (at declaration time)
     */
    Default(const std::string & value);

    /**
     * Factory function to make an empty default value which is obligatory.
     * This and following factory functions should be used instead of private constructors.
     *
     * Example of usage:
     * @code
     *      some_record.declare_key("some_key",Integer(),Default::obligatory(),"description");
     * @endcode
     */
    static Default obligatory()
    { return Default(obligatory_type); }

    /**
     * Factory function to make an empty default value which is optional.
     * To get the value of such key from the input you have to use non-throwing variant of the method
     * Input::Record::key, which returns the value through reference and allows checking presence of the key on the input.
     *
     * Example of usage:
     * @code
     *      some_record.declare_key("some_key",Integer(),Default::optional(),"description");
     * @endcode
     */
    static Default optional()
    { return Default(optional_type); }




    /**
     * Returns true if the default value should be specified at some time.
     * Currently, we support only default values given at declaration.
     */
    inline bool has_value() const
    { return (type_ == declaration); }

    inline bool is_obligatory() const
    { return (type_ == obligatory_type); }

    inline bool is_optional() const
    { return (type_ == optional_type); }

    /**
     * Returns stored value. Possibly empty string.
     */
    inline const string & value() const
    { return (value_); }

private:
    string value_;              ///< Stored value.
    enum DefaultType type_;     ///< Type of the Default.

    /**
     * Forbids default constructor.
     */
    Default();

    /**
     * Constructor for other types then 'declaration'.
     */
    Default(enum DefaultType type);

};



/**
 * @brief Base of classes for documentation and specification of input data.
 *
 * The main purpose of this class and its descendants is documentation and specification of the whole input tree.
 * It includes as the documentation of input structure for developers as documentation for users. Every class in the hierarchy
 * derived from @p TypeBase in the namespace @p Input::Type serves for description of types of input data (independently
 * on input file format). These types somehow follows C++ types: bool, int, double, string, enum, array, class since they are
 * intended to initialize these C++ data structures. In this analogy, declaration of descendants of @p TypeBase is like
 * specification of C++ language and instances of these types are like declarations of particular compound types -- classes
 * and arrays in C++ language.
 *
 * Basic scalar types
 *
 * The basic scalar types are String, Bool, Integer, Double, FileName. First four types directly corresponds to
 * C++ types. Individual instances of String and Bool are identical. On the other hand, instances of Integer and Double
 * can differ in interval of valid values. Finally, type FileName is  by @p filetype.
 *
 * Nontrivial scalar type is @p Selection<T>. It is template where @p T should be enum type that will be initialized from
 * data of this type. The Selection<T> object identify possible values of the enum type @p T with strings that should correspond to
 * names of the values. Unfortunately there is no way to get names of an enum type as strings so this information
 * has to be provided through method @p add_value.  Every @p Selection<T> object has particular name and description.
 * Description of individual values is optional.
 *
 * Array type
 *
 * Record type
 *
 * One instance of @p Input::Type::Record is like one class definition. You specify its members calling the method @p declare_key.
 * The keys should be valid C++ keywords. Every key represents a value of arbitrary type.
 *
 * Every class X that wants to initialize itself from the input should provide a static method  which returns
 * @p shared_ptr to the record used to initialize the class X. The shared pointer should be stored in a static
 * variable so that the record description is created only once for all instances of the class X.
 *
 * Types that are simple to initialize by calling their constructor (Array, Scalar) are cloned when used as a subtype in
 * Record::declare_key or Array::Array. On the other hand some other types have nontrivial initialization
 * (Record, AbstractRecord, Selection<T>). The latter group is not cloned on copy operations but rather use
 * shared_ptr to guarantee life time of the instance.
 *
 * AbstractRecord<class T>
 *
 * Mimics polymorphism of C++ classes. Data of this type can be used to initialize a variable of enum type T. The values of
 * this enum type are identified with particular Record types. These Records are like descendant of the AbstractRecord.
 *
 *  @ingroup input
 */
class TypeBase {
public:
    /**
     * Default constructor. Set all types finished after construction.
     */
    TypeBase() {}
    /**
     * @brief Implementation of documentation printing mechanism.
     *
     * It writes documentation into given stream @p stream. With @p extensive==false, it should print the description
     * about two lines long, while for @p extensive==true it outputs full documentation.
     * This is primarily used for documentation of Record types where we first
     * describe all keys of an Record with short descriptions and then we call recursively extensive documentation
     * for sub types that was not fully described yet. Further, we provide method @p reset_doc_flags to reset
     * all flags marking the already printed documentations. Parameter @p pad is used for correct indentation.
     *
     */
    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const = 0;

    /**
     * In order to output documentation of complex types only once, we mark types that have printed their documentation.
     * This method turns these marks off for the whole type subtree.
     */
    virtual void  reset_doc_flags() const =0;

    /// Returns true if the type is fully specified. In particular for Record and Selection, it returns true after @p finish() method is called.
    virtual bool is_finished() const
    {return true;}

    /// Returns an identification of the type. Useful for error messages.
    virtual string type_name() const  =0;

    virtual bool operator==(const TypeBase &other) const
        { return typeid(*this) == typeid(other); }

    bool operator!=(const TypeBase & other) const
        { return ! (*this == other); }

    /// Empty virtual destructor.
    virtual ~TypeBase( void ) {}

    virtual void valid_default(const string &str) const
    { }

protected:
    /**
     * Write out a string with given padding of every new line.
     */
    static std::ostream& write_description(std::ostream& stream, const string& str, unsigned int pad);

    /**
     * Type of hash values used in associative array that translates key names to indices in a record.
     *
     * For simplicity, we currently use whole strings as "hash".
     */
    typedef string KeyHash;

    /// Hash function.
    inline static KeyHash key_hash(const string &str) {
        return (str);
    }

    /**
     * Check that a @t key is valid identifier, i.e. consists only of valid characters, that are lower-case letters, digits and underscore,
     * we allow identifiers starting with a digit, but it is discouraged since it slows down parsing of the input file.
     */
    static bool is_valid_identifier(const string& key);


};

/**
 * For convenience we provide redirection operator for output documentation of Input:Type classes.
 */
std::ostream& operator<<(std::ostream& stream, const TypeBase& type);


class Record;
class SelectionBase;

/**
 * @ingroup input
 */
class Array : public TypeBase {

public:
    /**
     * Constructor with a @t type of array items given as pure reference. In this case @t type has to by descendant of @t TypeBase different from
     * 'complex' types @t Record and @ Selection<T>. You can also specify minimum and maximum size of the array.
     */
    template <class ValueType>
    inline Array(const ValueType &type, unsigned int min_size=0, unsigned int max_size=std::numeric_limits<unsigned int>::max() )
    : lower_bound_(min_size), upper_bound_(max_size)
    {
        // ASSERT MESSAGE: The type of declared keys has to be a class derived from TypeBase.
        BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, ValueType >::value) );

        ASSERT( min_size <= max_size, "Wrong limits for size of Input::Type::Array, min: %d, max: %d\n", min_size, max_size);

        if ( ! type.is_finished())
            xprintf(PrgErr, "Unfinished type '%s' used in declaration of Array.\n", type.type_name().c_str() );

        boost::shared_ptr<const ValueType> type_copy = boost::make_shared<ValueType>(type);
        type_of_values_ = type_copy;
        //set_type_impl(type, boost::is_base_of<TypeBase,ValueType>() );
    }


    /// Getter for the type of array items.
    inline const TypeBase &get_sub_type() const
        { return *type_of_values_; }

    /// Checks size of particular array.
    inline bool match_size(unsigned int size) const
        { return size >=lower_bound_ && size<=upper_bound_; }

    /// @brief Implements @p Type::TypeBase::documentation.
    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const;

    /// @brief Implements @p Type::TypeBase::reset_doc_flags.
    virtual void  reset_doc_flags() const;

    virtual string type_name() const;

    virtual bool operator==(const TypeBase &other) const
        { return  typeid(*this) == typeid(other) &&
                  (*type_of_values_ == static_cast<const Array *>(&other)->get_sub_type() );
        }

protected:



    boost::shared_ptr<const TypeBase> type_of_values_;
    unsigned int lower_bound_, upper_bound_;
};



/**
 * Base of scalar types.
 *
 *  @ingroup input
 */
class Scalar : public TypeBase {
public:

    virtual void  reset_doc_flags() const;
};


/**
 * @ingroup input
 */
class Bool : public Scalar {
public:
    Bool()
    {}

    virtual void valid_default(const string &str) const
    { from_default(str);}

    bool from_default(const string &str) const {
        if (str == "true" )  {
            return true;
        } else
        if (str == "false") {
            return false;
        } else {
            THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
        }
    }

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const;
    virtual string type_name() const;
};


/**
 * @ingroup input
 */
class Integer : public Scalar {
public:
    Integer(int lower_bound=std::numeric_limits<int>::min(), int upper_bound=std::numeric_limits<int>::max())
    : lower_bound_(lower_bound), upper_bound_(upper_bound)
    {}

    /**
     * Returns true if the given integer value conforms to the Type::Integer bounds.
     */
    bool match(int value) const {
        return ( value >=lower_bound_ && value <= upper_bound_);
    }

    /**
     * As before but also returns converted integer in @p value.
     */
    int from_default(const string &str) const {
        std::istringstream stream(str);
        int value;
        stream >> value;

        if (stream && stream.eof() && match(value)) {
            return value;
        } else {
            THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
        }
    }

    virtual void valid_default(const string &str) const
    { from_default(str);}



    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const;
    virtual string type_name() const;
private:
    int lower_bound_, upper_bound_;

};


/**
 * @ingroup input
 */
class Double : public Scalar {
public:
    Double(double lower_bound= -std::numeric_limits<double>::max(), double upper_bound=std::numeric_limits<double>::max())
    : lower_bound_(lower_bound), upper_bound_(upper_bound)
    {}

    /**
     * Returns true if the given integer value conforms to the Type::Double bounds.
     */
    bool match(double value) const {
        return ( value >=lower_bound_ && value <= upper_bound_);
    }


    /**
     * As before but also returns converted integer in @p value.
     */
    double from_default(const string &str) const {
        std::istringstream stream(str);
        double value;
        stream >> value;

        if (stream && stream.eof() && match(value)) {
            return value;
        } else {
            THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
        }
    }

    virtual void valid_default(const string &str) const
    { from_default(str);}


    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const;
    virtual string type_name() const;
private:
    double lower_bound_, upper_bound_;

};



/**
 * Just for consistency, but is essentialy same as Scalar.
 *
 * @ingroup input
 */
class String : public Scalar {
public:
    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const;
    virtual string type_name() const;

    virtual void valid_default(const string &str) const
    {from_default(str);}

    /**
     * Particular descendants can check validity of the string.
     */
    virtual string from_default(const string &str) const
    {return str;}
};


/**
 * @ingroup input
 */
class FileName : public String {
public:

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

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const;
    virtual string type_name() const;

    virtual bool operator==(const TypeBase &other) const
    { return  typeid(*this) == typeid(other) &&
                     (type_== static_cast<const FileName *>(&other)->get_file_type() );
    }

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
