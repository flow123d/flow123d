/*
 * InputType.hh
 *
 *  Created on: Mar 28, 2012
 *      Author: jb
 *
 *  todo:
 *
 *  - mozne rozdeleni na type.hh type_record_impl.hh type_selection_impl.hh type_impl.hh
 *
 *  - Doxygen doc
 *
 *  - projit vsechny tridy a zkontrolovat, zda poskytuji dostatecny pristup ke svemu obsahu.
 *
 *  - zavest dedeni Recordu (neco jako copy constructor), ale s tim, ze rodicovsky Record by si pamatoval sve potomky
 *    a tim bychom se zbavili AbstractRecordu -
 *    v kopii recordu by se mohly prepisovat klice (jak to zaridit, nemel by se z chyby udelat warning?
 *
 *  - zavest priznak pro record, ktery muze byt inicializovan z hodnoty jednoho klice (ten je treba zadat)
 *    a podobne pro pole.
 *    ...
 *
 *
 *  - ?? predelat DefaultValue na hierarchii trid (opet problem s error hlaskou) tkato by se to hlasilo pri kompilaci.
 *
 *  When C++11 specification become more supported, we can introduce class Key  that should be constructed form
 *  constant string during compilation, i particular it should check validity of the key string and compute the hash.
 *  This can provide some speedup for reading if it will be needed (probably not).
 */

#ifndef INPUTTYPE_HH_
#define INPUTTYPE_HH_

#include "system.hh"
#include "system/exceptions.hh"
#include <boost/type_traits.hpp>

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

/**
 * Macro to create a key object. The reason for this is twofold:
 * 1) We may use it to implement compile time computed hashes for faster access to the data.
 * 2) We can store line, function, filename where the key where used in order to report more specific error messages.
 */
#define KEY(name) #name

/**
 * Macro to simplify declaration of Record types. It declares local shared pointer @p var and initialize it
 * by new Record. The macro should be followed by parentheses with constructor parameters. Typical usage is:
 *
 * @code
 * MAKE_Input_Type_Record( transport_record )("TransportRecord", " Description of transport record.");
 * transport_record->declare_key(...);
 * @endcode
 */
#define MAKE_Input_Type_Record( var ) boost::shared_ptr<Record> var = boost::make_shared<Record>

/**
 * Macro to simplify declaration of Selection types. The first parameter is name of enum used as the template parameter the second
 * parameter of the macro is name of local shared pointer @p var to initialize by new Selection.
 * The macro should be followed by parentheses with constructor parameters. Typical usage is:
 *
 * @code
 * MAKE_Input_Type_Selection( enum Colors, colors_selection )("Colors");
 * colors_selection->add_value(...);
 * @endcode
 */
#define MAKE_Input_Type_Selection( enum_type, var) boost::shared_ptr< Selection<enum_type> > var = boost::make_shared< Selection<enum_type> >

namespace Input {
namespace Type {

using namespace std;


struct KeyNotFound : virtual FlowException {};
struct SelectionKeyNotFound : virtual FlowException {};
TYPEDEF_ERR_INFO( KeyName, const string );
TYPEDEF_ERR_INFO( RecordName, const string );
TYPEDEF_ERR_INFO( SelectionName, const string );

/**
 * @brief Possible file types.
 *
 *  TODO: move this declaration into possible global class Application
 *  and use also when opening the file and by IONAmeHandlar (which should be part of Application, and used by
 *  Input::Interface::Path (todo)
 *
 */
enum FileType {
    input_file,
    output_file
};


/**
 * @brief DefaultValue specifies default value of keys of a @p Record.
 *
 * It contains type of default value and possibly the value itself as a string.
 *
 * We prefer to use only default values specified at declaration since only those can by documented as part of
 * Record type specification.
 */
class DefaultValue {
public:
    /**
     * Possible types of default values.
     */
    enum DefaultType {
        declaration,    ///< Default value given at declaration time.
        optional,       ///< No default value, optional key. This is default type of the DefaultValue.
        obligatory      ///< No default value, obligatory key.
    };

    /**
     * Default constructor. Use type @t optional.
     */
    DefaultValue();

    /**
     * Constructor with given default value (at declaration time)
     */
    DefaultValue(const std::string & value);

    /**
     * Constructor for other types then 'declaration'.
     */
    DefaultValue(enum DefaultType type);

    /**
     * Returns true if the default value should be specified at some time.
     * Currently, we support only default values given at declaration.
     */
    inline bool has_value() const
    { return (type_ == declaration); }

    /**
     * Returns stored value. Possibly empty string.
     */
    inline const string & value() const
    { return (value_); }

private:
    string value_;              ///< Stored value.
    enum DefaultType type_;     ///< Type of the DefaultValue.
};


class Record;
class Scalar;
class Array;
class SelectionBase;


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
 */
class TypeBase {
public:
    /**
     * Default constructor. Set all types finished after construction.
     */
    TypeBase() : finished(true) {}
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
    inline bool is_finished() const
    {return finished;}

    /// Returns an identification of the type. Useful for error messages.
    virtual string type_name() const  =0;

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
    inline KeyHash key_hash(const string &str) const {
        return (str);
    }

    /**
     * Check that a @t key is valid identifier, i.e. consists only of valid characters, that are lower-case letters, digits and underscore,
     * we allow identifiers starting with a digit, but it is discouraged since it slows down parsing of the input file.
     */
    bool is_valid_identifier(const string& key);

    /// Empty virtual destructor.
    virtual ~TypeBase( void ) {}

    /**
     * Flag that is true for types with completed declaration. This provides nontrivial information only for Record and Selection.
     */
    bool finished;


};

/**
 * For convenience we provide redirection operator for output documentation of Input:Type classes.
 */
std::ostream& operator<<(std::ostream& stream, const TypeBase& type);


/**
 *
 * TODO: access to type_name_ ??
 */
class Record : public TypeBase {
private:
    /// Forbids usage of a copy constructor.
    Record(const Record& rec)
    {}

public:
    /**
     *  Structure for description of one key in record.
     *  The members dflt_type_ and default have reasonable meaning only for
     *  type_ == Scalar
     */
    struct Key {
        unsigned int key_index;                     ///< Position inside the record.
        string key_;                                ///< Key identifier.
        string description_;                        ///< Key description in context of particular Record type.
        boost::shared_ptr<const TypeBase> type_;    ///< Type of the key.
        DefaultValue default_;                      ///< DefaultValue, type and possibly value itself.
    };

    /**
     * Public typedef of constant iterator into array of keys.
     */
    typedef std::vector<struct Key>::const_iterator KeyIter;

    /**
     * Basic constructor. You has to provide @t type_name of the new declared Record type and
     * its @t description.
     */
    Record(const string & type_name, const string & description)
    : description_(description), type_name_(type_name), made_extensive_doc(false)
    {finished=false;}

    /**
     * Declares a key of the Record with name given by parameter @p key, the type given by parameter @p type, default value by parameter @p default_value, and with given
     * @p description. The parameter @p type has to be either boost::shared_ptr<TYPE> where TYPE is any of Record, AbstractRecord, and Selection, or @p type has to
     * be reference to any other descendant of TypeBase. The reason is that in the first group are types with complex description and would not to have multiple instances
     * of these types in the whole type hierarchy. This guarantees, that every Record, AbstractRecord, and Selection will be reported only once when documentation is printed out.
     */
    template <class KeyType>
    void declare_key(const string &key,
                            const KeyType &type,
                            const DefaultValue &default_value, const string &description)
    {
        if (finished) xprintf(PrgErr, "Declaration of key: %s in finished Record type: %s\n", key.c_str(), type_name_.c_str());

        if (! is_valid_identifier(key)) {
            xprintf(PrgErr, "Invalid key identifier %s in declaration of Record type: %s\n", key.c_str(), type_name_.c_str());
        }
        declare_key_impl(key,type, default_value, description, boost::is_base_of<TypeBase,KeyType>());
    }

    /**
     * Same as previous method but without given default value (same as DefaultValue(DefaultValue::none) )
     */
    template <class KeyType>
    void declare_key(const string &key,
                            const KeyType &type,
                            const string &description)
    {
        declare_key(key,type, DefaultValue(), description);
    }

    /**
     * Finish declaration of the Record type. Now further declarations can be added.
     */
    void finish() { finished = true; }

    /**
     * @brief Implements @p Type:TypeBase::documentation.
     */
    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const;

    /**
     * Set made_extensive_doc = false for this Record and all its descendants.
     */
    virtual void  reset_doc_flags() const;

    virtual string type_name() const;

    /**
     * Interface to mapping key -> index in record. Returns index (in continuous array) for given key.
     */
    inline unsigned int key_index(const string& key) const
    {
        ASSERT( finished, "Asking for information of unfinished Record type: %s\n", type_name_.c_str());
        KeyHash key_h = key_hash(key);
        key_to_index_const_iter it = key_to_index.find(key_h);
        if (it != key_to_index.end()) return it->second;
        else
            throw KeyNotFound() << KeyName_EI(key) << RecordName_EI(type_name_);

        return size();
    }

    /**
     * Returns iterator to the key struct for given key string.
     *
     */
    inline KeyIter key_iterator(const string& key) const
    {
        ASSERT( finished, "Asking for information of unfinished Record type: %s\n", type_name_.c_str());
        return begin() + key_index(key);
    }

    inline KeyIter begin() const
    {
        ASSERT( finished, "Asking for information of unfinished Record type: %s\n", type_name_.c_str());
        return keys.begin();
    }

    inline KeyIter end() const
    {
        ASSERT( finished, "Asking for information of unfinished Record type: %s\n", type_name_.c_str());
        return keys.end();
    }

    inline int size() const {
        ASSERT( finished, "Asking for information of unfinished Record type: %s\n", type_name_.c_str());
        ASSERT( keys.size() == key_to_index.size(), "Sizes of Type:Record doesn't match. (map: %d vec: %d)\n", key_to_index.size(), keys.size());
        return keys.size();
    }


protected:
    /**
     * Ultimate implementation of all decalare_key methods. Using just boost::shared_ptr for @p type.
     */
    template <class KeyType>
    void declare_key_impl(const string &key,
                     boost::shared_ptr<KeyType> type,
                     const DefaultValue &default_value, const string &description,
                     const boost::false_type &)
    {
        // ASSERT MESSAGE: The type of declared keys has to be a class derived from TypeBase.
        BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, KeyType>::value) );

        if ( (! type->is_finished()) &&  ( (void*)(type.get()) != (void*)(this) ) )
            xprintf(PrgErr, "Unfinished type of declaring key: %s in Record type: %s\n", key.c_str(), type_name_.c_str() );

        // If KeyType is not derived from Scalar, we check emptiness of the default value.
        if (boost::is_base_of<Scalar, KeyType>::value == false && default_value.has_value() ) {
            xprintf(Err, "Default value for non scalar type in declaration of key: %s in Record type: %s \n", key.c_str(), type_name_.c_str() );
        }

        KeyHash key_h = key_hash(key);
        if (key_to_index.find(key_h) == key_to_index.end()) {
           key_to_index.insert( std::make_pair(key_h, keys.size()) );
           // make our own copy of type object allocated at heap (could be expensive, but we don't care)
           Key tmp_key = { (unsigned int)keys.size(), key, description, type, default_value};
           keys.push_back(tmp_key);
        } else {
           xprintf(Err,"Re-declaration of the key: %s in Record type: %s\n", key.c_str(), type_name_.c_str() );
        }
    }

    /**
     * Implementation for types given by reference - creates auxiliary shared_ptr.
     */
    template <class KeyType>
    inline void declare_key_impl(const string &key,
            const KeyType &type,
            const DefaultValue &default_value, const string &description, const boost::true_type &)
    {
        /// ASSERT MESSAGE: "You have to use shared_ptr to declare key with types Record or Selection. For those classes we forbids copies."
        BOOST_STATIC_ASSERT( ( ! boost::is_same<Record,KeyType>::value &&  ! boost::is_base_of<SelectionBase,KeyType>::value) );
        boost::shared_ptr<const KeyType> type_copy = boost::make_shared<KeyType>(type);
        declare_key_impl(key,type_copy, default_value, description, boost::false_type());
    }




protected:
    /// Database of valid keys
    std::map<KeyHash, unsigned int> key_to_index;
    typedef std::map<KeyHash, unsigned int>::const_iterator key_to_index_const_iter;

    /// Keys in order as they where declared.
    std::vector<struct Key> keys;


    /// Description of the whole record type.
    const string description_;
    const string type_name_;

    /**
     * This flag is set to true when documentation of the Record was called with extensive==true
     * and full description of the Record was produced.
     *
     * This member is marked 'mutable' since it doesn't change structure or description of the type. It only influence the output.
     */
    mutable bool made_extensive_doc;

};





/**
 *
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
        ASSERT( min_size <= max_size, "Wrong limits for size of Input::Type::Array, min: %d, max: %d\n", min_size, max_size);
        // ASSERT MESSAGE: The type of declared keys has to be a class derived from TypeBase.
        BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, ValueType >::value) );
        /// ASSERT MESSAGE: "You have to use shared_ptr to declare key with types Record or Selection. For those classes we forbids copies."
        BOOST_STATIC_ASSERT( ( ! boost::is_same<Record,ValueType>::value &&  ! boost::is_base_of<SelectionBase, ValueType >::value) );

        if ( ! type.is_finished())
            xprintf(PrgErr, "Unfinished type '%s' used in declaration of Array.\n", type.type_name().c_str() );

        boost::shared_ptr<const ValueType> type_copy = boost::make_shared<ValueType>(type);
        type_of_values_ = type_copy;
        //set_type_impl(type, boost::is_base_of<TypeBase,ValueType>() );
    }

    /**
     * Constructor with a @t type of array items given through shared_ptr. In this case @t type has to by descendant of @t TypeBase.
     * You can also specify minimum and maximum size of the array.
     */
    template <class ValueType>
    Array(boost::shared_ptr<ValueType> type, unsigned int min_size=0, unsigned int max_size=std::numeric_limits<unsigned int>::max() )
    : type_of_values_(type),lower_bound_(min_size), upper_bound_(max_size)
    {
        if ( ! type->is_finished())
            xprintf(PrgErr, "Unfinished type '%s' used in declaration of Array.\n", type->type_name().c_str() );

        // ASSERT MESSAGE: The type of declared keys has to be a class derived from TypeBase.
        BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, ValueType>::value) );
    }



    /// Getter for the type of array items.
    inline const TypeBase &get_sub_type() const
        { return *type_of_values_; }

    /// @brief Implements @p Type::TypeBase::documentation.
    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const;

    /// @brief Implements @p Type::TypeBase::reset_doc_flags.
    virtual void  reset_doc_flags() const;

    virtual string type_name() const;
protected:



    boost::shared_ptr<const TypeBase> type_of_values_;
    unsigned int lower_bound_, upper_bound_;
};



/**
 * Base of scalar types.
 */
class Scalar : public TypeBase {
public:

    virtual void  reset_doc_flags() const;
};

/**
 * Base of Selection templates.
 */
class SelectionBase : public Scalar {
};

/**
 * In order to by type safe this class has to be templated by actual enum to which
 * the input value could be converted.
 *
 * This template assumes that Enum is an enum type or enum class type (in C++11).
 *
 * It would be nice to have specialization that could be constructed from enum types produced by CppEnumMacro.
 * Then we can drop add_value method.
 *
 * Future shows, if this is not too restrictive. Maybe, it is more practical drop the template and use just plain ints for values.
 */
template <class Enum>
class Selection : public SelectionBase {
public:
    Selection(const string &name)
    :type_name_(name), made_extensive_doc(false)
    { finished = false;}

    /**
     * Adds one new @p value with name given by @p key to the Selection. The @p description of meaning of the value could be provided.
     */
    void add_value(const Enum value, const std::string &key, const std::string &description = "") {
        F_ENTRY;

        if (finished) xprintf(PrgErr, "Declaration of new name: %s in finished Selection type: %s\n", key.c_str(), type_name_.c_str());
        KeyHash key_h = key_hash(key);
        if (key_to_index_.find(key_h) != key_to_index_.end()) {
            xprintf(PrgErr,"Name '%s' already exists in Selection: %s\n", key.c_str(), type_name_.c_str());
            return;
        }
        value_to_index_const_iter it = value_to_index_.find(value);
        if ( it  != value_to_index_.end()) {
            xprintf(PrgErr,"Value %d of new name '%s' conflicts with value %d of previous name '%s' in Selection: '%s'.\n",
                     value, key.c_str(), keys_[it->second].value, keys_[it->second].key_.c_str(), type_name_.c_str());
            return;
        }

        unsigned int new_idx= key_to_index_.size();
        key_to_index_.insert( std::make_pair(key_h, new_idx) );
        value_to_index_.insert( std::make_pair(value, new_idx) );

        Key tmp_key = { new_idx, key, description, value};
        keys_.push_back(tmp_key);
    }

    void finish() {finished=true;}

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const {
        if (!finished)  xprintf(PrgErr, "Can not provide documentation of unfinished Selection type: %s\n", type_name_.c_str());

        if (!extensive) {
            stream << "Selection '"<< type_name_ << "' of " << size() << " values.";
        }
        if (extensive && !made_extensive_doc) {
            made_extensive_doc=true;

            stream << endl << "Selection '"<< type_name_ << "' of " << size() << " values." << endl;
            stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
            // keys
            for(keys_const_iterator it = keys_.begin(); it!=keys_.end(); ++it) {
                stream << setw(4) << ""
                       << it->key_ << " = " << it->value;
                if (it->description_ != "") stream  << " (" << it->description_ << ")";
                stream << endl;
            }
            stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ')
                        << " " << type_name_ << endl;
        }
        return stream;
    }


    virtual void  reset_doc_flags() const {
        made_extensive_doc=false;
    }

    virtual string type_name() const {
        return type_name_;
    }

    /***
     * Converts name (on input) to the value. Throws if the name do not exist.
     */
    inline Enum name_to_value(const string &key)
    {
        ASSERT( finished, "Asking for information of unfinished Selection type: %s\n", type_name_.c_str());
        KeyHash key_h = key_hash(key);
        key_to_index_const_iter it = key_to_index_.find(key_h);
        if (it != key_to_index_.end()) return ( keys_[it->second].value );
        else
            throw SelectionKeyNotFound() << KeyName_EI(key) << SelectionName_EI(type_name_);

     }

    /**
     * Just check if there is a particular name in the Selection.
     */
    inline bool has_name(const string &key) {
        ASSERT( finished, "Asking for information of unfinished Selection type: %s\n", type_name_.c_str());
        KeyHash key_h = key_hash(key);
        return  ( key_to_index_.find(key_h) != key_to_index_.end() ) ;
    }

    /***
     *  Check if there is a particular value in the Selection.
     */
    inline bool has_value(Enum val) {
        ASSERT( finished, "Asking for information of unfinished Selection type: %s\n", type_name_.c_str());
        return ( value_to_index_.find(val) != value_to_index_.end() );
    }


    inline unsigned int size() const {
        ASSERT( finished, "Asking for information of unfinished Selection type: %s\n", type_name_.c_str());
        ASSERT( keys_.size() == key_to_index_.size(), "Sizes of Type:Selection doesn't match. (map: %d vec: %d)\n", key_to_index_.size(), keys_.size());
        return keys_.size();
    }

private:
    string type_name_;

    struct Key {
        unsigned int key_index;
        string key_;
        string description_;
        Enum value;
    };
    /// Map of valid keys to index.
    std::map<KeyHash, unsigned int> key_to_index_;
    typedef std::map<KeyHash, unsigned int>::const_iterator key_to_index_const_iter;

    /// Map of valid values to index.
    typename std::map<Enum, unsigned int> value_to_index_;
    typedef typename std::map<Enum, unsigned int>::const_iterator value_to_index_const_iter;

    std::vector<Key> keys_;
    typedef typename std::vector<struct Key>::const_iterator keys_const_iterator;

    /**
     * This flag is set to true when documentation of the Record was called with extensive==true
     * and full description of the Record was produced.
     *
     * This member is marked 'mutable' since it doesn't change structure or description of the type. It only influence the output.
     */
    mutable bool made_extensive_doc;

};

/**
 *
 */
class Bool : public Scalar {
public:
    Bool()
    {}

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const;
    virtual string type_name() const;
};


/**
 *
 */
class Integer : public Scalar {
public:
    Integer(int lower_bound=std::numeric_limits<int>::min(), int upper_bound=std::numeric_limits<int>::max())
    : lower_bound_(lower_bound), upper_bound_(upper_bound)
    {}

    bool match(string &str) {
        int value;
        return match(str, value);
    }

    bool match(string &str, int &value) {
        std::istringstream stream(str);
        stream >> value;
        if (stream.good()) {
            return match(value);
        } else {
            return false;
        }
    }

    bool match(int value) {
        return ( value >=lower_bound_ && value <= upper_bound_);
    }

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const;
    virtual string type_name() const;
private:
    int lower_bound_, upper_bound_;

};


/**
 *
 */
class Double : public Scalar {
public:
    Double(double lower_bound=std::numeric_limits<double>::min(), double upper_bound=std::numeric_limits<double>::max())
    : lower_bound_(lower_bound), upper_bound_(upper_bound)
    {}


    bool match(string &str) {
        double value;
        return match(str, value);
    }

    bool match(std::string &str, double &value) {
        std::istringstream stream(str);
        stream >> value;
        if (stream.good()) {
            return match(value);
        } else {
            return false;
        }
    }

    bool match(double value) {
        return ( value >=lower_bound_ && value <= upper_bound_);
    }

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const;
    virtual string type_name() const;
private:
    double lower_bound_, upper_bound_;

};



/**
 * Just for consistency, but is essentialy same as Scalar.
 */
class String : public Scalar {
public:
    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const;
    virtual string type_name() const;
};


/**
 *
 */
class FileName : public String {
public:
    FileName(FileType type)
    : type_(type)
    {}

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const;
    virtual string type_name() const;

    FileType get_file_type() const {
        return type_;
    }



private:
    FileType    type_;
};

} // closing namespace Type
} // closing namespace Input


#endif /* INPUTTYPE_HH_ */
