/*
 * InputType.hh
 *
 *  Created on: Mar 28, 2012
 *      Author: jb
 *
 *  todo:
 *
 *  - Check lower case and no whitespace in key strings.
 *  - zakazat copy constructor v Record - testovat, ze nejde volat declare_key s Recordem mimo shared_ptr
 *  - zavest dedeni Recordu (neco jako copy constructor), ale s tim, ze rodicovsky Record by si pamatoval sve potomky
 *    a tim bychom se zbavili AbstractRecordu -
 *    v kopii recordu by se mohly prepisovat klice (jak to zaridit, nemel by se z chyby udelat warning?
 *  - zavest uzavreni Recordu a Selection
 *  - zavest priznak pro record, ktery muze byt inicializovan z hodnoty jednoho klice (ten je treba zadat)
 *    a podobne pro pole.
 *    ...
 *
 *  - dokumentace hotovych trid
 *  - presun neinline metod do *.cc
 *
 *  - ? predelat DefaultValue na hierarchii trid (opet problem s error hlaskou) tkato by se to hlasilo pri kompilaci.
 */

#ifndef INPUTTYPE_HH_
#define INPUTTYPE_HH_

#include <limits>
#include <ios>
#include <map>
#include <vector>
#include <string>
#include <iomanip>

#include "system.hh"

#include <boost/type_traits.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

/**
 * Macro to create a key object. The reason for this is twofold:
 * 1) We may use it to implement compile time computed hashes for faster access to the data.
 * 2) We can store line, function, filename where the key where used to
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

// TODO: move FileType into Application and use also when opening the file and by IONAmeHandlar
// (which in turn should be part of input interface)
enum FileType {
    input_file,
    output_file
};


/**
 * @brief DefaultValue specifies type and possibly string of default value used for keys of a @p Record.
 */
class DefaultValue {
public:
    /**
     * Possible types of default values.
     */
    enum DefaultType {
        none,           ///< no specification of default value
        read_time,      ///< default value will be given at read time
        declaration,    ///< default value given at declaration time (can be overwritten at read time)
        optional,       ///< no default value, optional key
        obligatory      ///< no default value, obligatory key
    };

    /**
     * Default constructor provides NONE default value.
     */
    DefaultValue()
    : value_(), type_(none)
    {}

    /**
     * Constructor with given default value (at declaration time)
     */
    DefaultValue(const std::string & value)
    : value_(value), type_(declaration)
    {}

    /**
     * Constructor for other types then 'declaration'.
     */
    DefaultValue(enum DefaultType type)
    : value_(), type_(type)
    {
        if (type == declaration) {
            xprintf(Err, "Can not construct DefaultValue with type 'declaration' without providing the default value.\n");
        }
    }

    /**
     * Returns true if the default value should be specified at some time.
     */
    bool given_value() const
    { return (type_ == declaration || type_ == read_time); }

    /**
     * Returns stored value. Possibly empty string.
     */
    const string & value() const
    { return (value_); }

private:
    string value_;
    enum DefaultType type_;
};

class AbstractRecord;
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
     * Fundamental method for output documentation of Input Types.
     * It writes documentation into given stream @param stream. Parameter @param deep
     * indicates verbosity of produced documentation. With extensive==false, the description should
     * be about two lines. This is primarily used for documentation of Record types where we first
     * describe all keys of an Record with short descriptions and then call recursively extensive documentation
     * for sub types that was not fully described yet. Further we provide static function to reset
     * all flags marking state of documentation.
     *
     */
    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const = 0;

    /**
     * In order to output documentation of complex types only once, we mark types that have printed their documentation.
     * This method turns these marks off for the whole type subtree.
     */
    virtual void  reset_doc_flags() const {
    }

    /**
     * Write out a string with given padding of every new line.
     */
    static std::ostream& write_description(std::ostream& stream, const string& str, unsigned int pad) {
        boost::tokenizer<boost::char_separator<char> > line_tokenizer(str, boost::char_separator<char>("\n"));
        boost::tokenizer<boost::char_separator<char> >::iterator tok;

        // Up to first \n without padding.
        stream << endl;

            // For every \n add padding at beginning of the nex line.
            for(tok = line_tokenizer.begin(); tok != line_tokenizer.end(); ++tok) {
                stream << setw(pad) << "" << "# "
                        << *tok << endl;
            }
        return stream;
    }

    /**
     * For simplicity we use key strings directly as keys in the associative array (map or hash table).
     * However hashes may be used to speedup the lookup. To this end we declare here type
     * of keys in associative array and a "hashing function" to produce these values from strings.
     */
    typedef string KeyHash;

    /// Envelop around compile time hashing macro. To provide same hashing at runtime.
    KeyHash key_hash(const string &str) const {
        return (str);
    }

    virtual ~TypeBase( void ) {}
};

/**
 * For convenience we provide redirection operator for output documentation of Input:Type classes.
 */
std::ostream& operator<<(std::ostream& stream, const TypeBase& type);

class AbstractRecord : public TypeBase {

};


/**
 *
 *
 */
class Record : public TypeBase {
public:
    /**
     *  Structure for description of one key in record.
     *  The members dflt_type_ and default have reasonable meaning only for
     *  type_ == Scalar
     */
    ///
    struct Key {
        unsigned int key_index;
        string key_;
        string description_;
        boost::shared_ptr<const TypeBase> type_;
        DefaultValue default_;
    };

    /**
     * Public typedef of constant iterator into array of keys.
     */
    typedef std::vector<struct Key>::const_iterator KeyIter;

    Record(const string & type_name, const string & description)
    : description_(description), type_name_(type_name), made_extensive_doc(false)
    {}

    /**
     * Declares a key of the Record with name given by parameter @p key, the type given by parameter @p type, default value by parameter @p default_value, and with given
     * @p description. The parameter @p type has to be either boost::shared_ptr<TYPE> where TYPE is any of Record, AbstractRecord, and Selection, or @p type has to
     * be reference to any other descendant of TypeBase. The reason is that in the first group are types with complex description and would not to have multiple instances
     * of these types in the whole type hierarchy. This guarantees, that every Record, AbstractRecord, and Selection will be reported only once when documentation is printed out.
     */
    template <class KeyType>
    inline void declare_key(const string &key,
                            const KeyType &type,
                            const DefaultValue &default_value, const string &description)
    {
        declare_key_impl(key,type, default_value, description, boost::is_base_of<TypeBase,KeyType>());
    }

    /**
     * Same as previous method but without given default value (same as DefaultValue(DefaultValue::none) )
     */
    template <class KeyType>
    inline void declare_key(const string &key,
                            const KeyType &type,
                            const string &description)
    {
        declare_key(key,type, DefaultValue(), description);
    }

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const {

        if (! extensive) {

            // Short description
            stream << "Record '" << type_name_ << "' with "<< keys.size() << " keys";
        } else if ( ! made_extensive_doc) {

            // Extensive description
            made_extensive_doc=true;

            // header
            stream << endl;
            stream << ""
                   << "Record '" << type_name_ << "' with "<< keys.size() << " keys.";
            write_description(stream, description_, pad);
            stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
            // keys
            for(KeyIter it = keys.begin(); it!=keys.end(); ++it) {
                stream << setw(pad + 4) << ""
                       << it->key_ << " = <" << it->default_.value() << "> is ";
                it->type_->documentation( stream , false, pad +4 ); // short description of the type of the key
                write_description(stream, it->description_, pad + 4); // description of the key on further lines

            }
            stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ')
            << " " << type_name_ << endl;

            // Full documentation of embedded record types.
            for(KeyIter it = keys.begin(); it!=keys.end(); ++it) {
                it->type_->documentation(stream, true, 0);
            }

        }

        return stream;
    }



    /**
     * Set made_extensive_doc = false for this Record and all its descendants.
     */
    virtual void  reset_doc_flags() {
        made_extensive_doc=false;
        for(KeyIter it = keys.begin(); it!=keys.end(); ++it) {
            it->type_->reset_doc_flags();
        }
    }

    /**
     * Interface to mapping key -> index in record. Returns index (in continuous array) for given key.
     *
     * TODO: Throw exception instead of message, to report more precise error message at higher level.
     *
     */
    inline unsigned int key_index(const string& key) const
    {
        KeyHash key_h = key_hash(key);
        key_to_index_const_iter it = key_to_index.find(key_h);
        if (it != key_to_index.end()) return it->second;
        else
            xprintf(Err, "Attempt to read key '%s', which is not declared within Record '%s'\n", key.c_str(), type_name_.c_str() );

        return size();
    }

    /**
     * Returns iterator to the key struct for given key string.
     *
     */
    KeyIter key_iterator(const string& key) const
    {
        return begin() + key_index(key);
    }

    inline KeyIter begin() const
    {
        return keys.begin();
    }

    inline KeyIter end() const
    {
        return keys.end();
    }

    inline int size() const {
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
        // catch call with shared_ptr to some other type
        BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, KeyType>::value) );
        // check if KeyType is derived from Scalar so that default value can be given.
        if (boost::is_base_of<Scalar, KeyType>::value == false && default_value.given_value() ) {
            xprintf(Err, "Can not provide default value for non scalar type. Key: %s\n", key.c_str());
        }

        KeyHash key_h = key_hash(key);
        if (key_to_index.find(key_h) == key_to_index.end()) {
           key_to_index.insert( std::make_pair(key_h, keys.size()) );
           // make our own copy of type object allocated at heap (could be expensive, but we don't care)
           Key tmp_key = { (unsigned int)keys.size(), key, description, type, default_value};
           keys.push_back(tmp_key);
        } else {
           xprintf(Err,"Redeclaration of key: %s\n", key.c_str());
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
        if (boost::is_same<Record,KeyType>() || boost::is_base_of<SelectionBase,KeyType>()) {
            xprintf(Err,"Complex type (Record, Selection, AbstractRecord) has to by passed as shared_ptr. Declaration of key: %s\n", key.c_str());
        }
        boost::shared_ptr<const KeyType> type_copy = boost::make_shared<KeyType>(type);
        declare_key_impl(key,type_copy, default_value, description, boost::false_type());
    }

private:

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
private:
    boost::shared_ptr<const TypeBase> type_of_values_;
    unsigned int lower_bound_, upper_bound_;

    template <class ValueType>
    void set_type_impl(boost::shared_ptr<ValueType> type,  const boost::false_type &)
    {
        type_of_values_ = type;
    }

    template <class ValueType>
    inline void set_type_impl(const ValueType &type,  const boost::true_type &)
    {
        if (boost::is_same<Record,ValueType>()) {
            xprintf(Err,"Complex type (Record, Selection, AbstractRecord) has to by passed as shared_ptr.\n");
        }
        boost::shared_ptr<const ValueType> type_copy = boost::make_shared<ValueType>(type);
        set_type_impl(type_copy, boost::false_type());
    }


public:
    template <class ValueType>
    inline Array(const ValueType &type, unsigned int min_size=0, unsigned int max_size=std::numeric_limits<unsigned int>::max() )
    : lower_bound_(min_size), upper_bound_(max_size)
    {
        set_type_impl(type, boost::is_base_of<TypeBase,ValueType>() );
    }

    template <class ValueType>
    Array(boost::shared_ptr<const ValueType> type, unsigned int min_size=0, unsigned int max_size=std::numeric_limits<unsigned int>::max() )
    : type_of_values_(type),lower_bound_(min_size), upper_bound_(max_size)
    {}

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const {
        if (extensive) {
            type_of_values_->documentation(stream, true, pad);
        } else {
            stream << "Array, size limits: [" << lower_bound_ << ", " << upper_bound_ << "] of type: " << endl;
            stream << setw(pad+4) << "";
            type_of_values_->documentation(stream, false, pad+4);
        }
        return stream;
    }

    virtual void  reset_doc_flags() const {
        type_of_values_->reset_doc_flags();
    }

    const TypeBase &get_sub_type() const {
        return *type_of_values_;
    }

};



/**
 * Base of scalar types.
 */
class Scalar : public TypeBase {
public:

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const {

        if (extensive) return stream;
        stream << "String (generic)";
        return stream;
    }
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
 * It would be nice to have specialization that could be constructed from enum types procuced by CppEnumMacro.
 * Then we can drop add_value method.
 */
template <class Enum>
class Selection : public SelectionBase {
public:
    Selection(const string &name)
    :name_(name), made_extensive_doc(false)
    {}

    void add_value(const Enum value, const std::string &key, const std::string &description = "") {
        F_ENTRY;

        KeyHash key_h = key_hash(key);
        if (key_to_index_.find(key_h) != key_to_index_.end()) {
            xprintf(Err,"Existing name in declaration of enum key, name: %s value: %d\n", key.c_str(), value);
            return;
        }
        if (value_to_index_.find(value) != value_to_index_.end()) {
            xprintf(Err,"Existing value in declaration of enum key, name: %s value: %d\n", key.c_str(), value);
            return;
        }

        unsigned int new_idx= size();
        key_to_index_.insert( std::make_pair(key_h, new_idx) );
        value_to_index_.insert( std::make_pair(value, new_idx) );

        Key tmp_key = { new_idx, key, description, value};
        keys_.push_back(tmp_key);
    }

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const {
        if (!extensive) {
            stream << "Selection '"<< name_ << "' of " << size() << " values.";
        }
        if (extensive && !made_extensive_doc) {
            made_extensive_doc=true;

            stream << endl << "Selection '"<< name_ << "' of " << size() << " values." << endl;
            stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
            // keys
            for(keys_const_iterator it = keys_.begin(); it!=keys_.end(); ++it) {
                stream << setw(4) << ""
                       << it->key_ << " = " << it->value;
                if (it->description_ != "") stream  << " (" << it->description_ << ")";
                stream << endl;
            }
            stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ')
                        << " " << name_ << endl;
        }
        return stream;
    }

    inline unsigned int size() const {
        ASSERT( keys_.size() == key_to_index_.size(), "Sizes of Type:Selection doesn't match. (map: %d vec: %d)\n", key_to_index_.size(), keys_.size());
        return keys_.size();
    }

    virtual void  reset_doc_flags() const {
        made_extensive_doc=false;
    }

private:
    string name_;

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

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const {
        if (extensive) return stream;
        stream << "Bool";
        return stream;
    }
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

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const {
        if (extensive) return stream;
        stream << "Integer in [" << lower_bound_ << ", " << upper_bound_ << "]";
        return stream;
    }
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

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const {
        if (extensive) return stream;

        stream << "Double in [" << lower_bound_ << ", " << upper_bound_ << "]";
        return stream;
    }
private:
    double lower_bound_, upper_bound_;

};



/**
 * Just for consistency, but is essentialy same as Scalar.
 */
class String : public Scalar {
};


/**
 *
 */
class FileName : public String {
public:
    FileName(FileType type)
    : type_(type)
    {}

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const {
        if (extensive) return stream;

        stream << "FileName of ";
        switch (type_) {
        case input_file:
            stream << "input file";
            break;
        case output_file:
            stream << "output file";
            break;
        default:
            stream << "file with unknown type";
            break;
        }
        return stream;
    }

private:
    FileType    type_;
};

} // closing namespace Type
} // closing namespace Input


#endif /* INPUTTYPE_HH_ */
