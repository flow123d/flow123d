/*
 * InputType.hh
 *
 *  Created on: Mar 28, 2012
 *      Author: jb
 *
 *  todo:
 *
 *  Check lower case and no whitespace in key strings.
 *
 *  Jak zaridit, aby se v Recodr. declare_key a v konstruktoru Array, vytvareli kopie jen
 *  pro typy != Record.
 *
 *
 *  -# budu mit virtulani metodu clone(), ktera pro Record, AbstractRecord a Enum bude vracet kopii jejich shared pointeru
 *     pro ostatni provede skutecnou kopii
 *  -# tj. i vsechno predavani typu musi byt pres shared_ptr
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

namespace Input {
namespace Type {


using std::string;

// TODO: move FileType into Application and use also when opening the file and by IONAmeHandlar
// (which in turn should be part of input interface)
enum FileType {
    input_file,
    output_file
};

class DefaultValue {
public:
    /// This enum says when the default value shoud be specified.
    enum DefaultType {
        none,
        read_time,
        declaration
    };

    /**
     * Default constructor provides NONE default value.
     */
    DefaultValue()
    : value_(), type_(none)
    {}

    /**
     * Constructor with given default value.
     */
    DefaultValue(const std::string & value)
    : value_(value), type_(declaration)
    {}
    /**
     * Constructor for default given at read time.
     */
    DefaultValue(enum DefaultType type)
    : value_(), type_(type)
    {
        if (type == declaration) {
            xprintf(Err, "Can not construct DefaultValue with type 'declaration' without providing the default value.\n");
        }
    }


    bool not_none() const
    { return (type_ != none); }

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


/**
 * @brief Base of classes for documentation and specification of input data.
 *
 * There are three groups of input data:
 * -# scalar data represented by an unstructured  string
 * -# list of ordered data of a same type, can be indexed
 * -# record of data of various type, indexed by string keys
 *
 * A record type should be declared in the static method of the class that reads from this record. This
 * static method should have a static variable -- pointer to the Record type created at first call and
 * just returned latter on.
 *
 * The Record type stores copy of type of every individual key. For very big inputs one can consider using boost:flyweight
 * in order to reduce redundancy.
 *
 * Types that are simple to initialize by calling their constructor (Array, Scalar) are cloned when used as a subtype in
 * Record::declare_key or Array::Array. On the other hand some other types are named with nontrivial initialization (Record, AbstractRecord, Selection).
 * The latter group is not cloned on copy operations but rather use shared_ptr to guarantee life time of the instance.
 *
 * For the first group we assume usage of temporary objects while for the latter group all  construcors are private, but we provide static factory methods
 * to guarantee that all abjects are referenced through shared_ptr.
 */
class TypeBase {
protected:


public:
    /**
     * Fundamental method for output documentation of Input Types.
     * It writes documentation into given stream @param stream. Parameter @param deep
     * indicates verbosity of produced documentation. With extensive==false, the description should
     * be about two lines. This is primarly used for documentation of Record types where we first
     * describe all keys of an Record with short descriptions and then call recursively extensive documentation
     * for sub types that was not fully described yet. Further we provide static function to reset
     * all flags marking state of documentation.
     *
     */
    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const = 0;
    /**
     *
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
};

/**
 * For convenience we provide redirection operator for output documentation of Input:Type classes.
 */
std::ostream& operator<<(std::ostream& stream, const TypeBase& type) {
    return type.documentation(stream);
}

class AbstractRecord : public TypeBase {

};



/**
 * Possible types of default values:
 * - given at declaration time
 * - provided at read time (not part of Record specification)
 * - none
 *
 *
 */
class Record : public TypeBase {

public:
    //static inline boost::shared_ptr<Record> & record_factory(const string & type_name, const string & description)
    //{
    //    return boost::make_shared(Record(type_name, description) );
    //}
    Record(const string & type_name, const string & description)
    : type_name_(type_name), description_(description), made_extensive_doc(false)
    {}


protected:
    /**
     * Method for declaration of keys of scalar types with given default value.
     *
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
        if (boost::is_base_of<Scalar, KeyType>::value == false && default_value.not_none() ) {
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

    template <class KeyType>
    inline void declare_key_impl(const string &key,
            const KeyType &type,
            const DefaultValue &default_value, const string &description, const boost::true_type &)
    {
        boost::shared_ptr<const KeyType> type_copy = boost::make_shared<KeyType>(type);
        declare_key_impl(key,type_copy, default_value, description, boost::false_type());
    }

public:
    /**
     * Version without given default value.
     */
    template <class KeyType>
    inline void declare_key(const string &key,
                            const KeyType &type,
                            const DefaultValue &default_value, const string &description)
    {
        declare_key_impl(key,type, default_value, description, boost::is_base_of<TypeBase,KeyType>());
    }

    /**
     * Version without given default value.
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
            stream << "Record '" << type_name_ << "' with "<< keys.size() << " keys" << std::endl;
        } else if ( ! made_extensive_doc) {

            // Extensive description
            made_extensive_doc=true;
            // header
            stream << endl;
            stream << setw(pad) << ""
                   << "Record '" << type_name_ << "' with "<< keys.size() << " keys.";
            write_description(stream, description_, pad) << endl;

            // keys
            for(keys_const_iterator it = keys.begin(); it!=keys.end(); ++it) {
                stream << setw(pad + 4) << ""
                       << it->key_ << " = <" << it->default_.value() << "> is ";
                it->type_->documentation( stream , false, pad +4 ); // short description of the type of the key
                write_description(stream, it->description_, pad + 4); // description of the key on further lines
                stream << endl;
            }

            // Full documentation of embedded record types.
            for(keys_const_iterator it = keys.begin(); it!=keys.end(); ++it) {
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
        for(keys_const_iterator it = keys.begin(); it!=keys.end(); ++it) {
            it->type_->reset_doc_flags();
        }
    }

    //virtual const &TypeBase clone() const {
    //    return
    //}


private:


    /**
     * For simplicity we use key strings directly as keys in the associative array (map or hash table).
     * However hashes may be used to speedup the lookup. To this end we declare here type
     * of keys in associative array and a "hashing function" to produce these values from strings.
     */
    typedef string KeyHash;
    /// Envelop around compile time hashing macro. To provide same hashing at runtime.
    KeyHash key_hash(const string &str) {
        return (str);
    }


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

    /// Database of valid keys.
    std::map<KeyHash, unsigned int> key_to_index;

    /// Keys in order as they where declared.
    std::vector<struct Key> keys;
    typedef std::vector<struct Key>::const_iterator keys_const_iterator;

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


class Array : public TypeBase {
private:
    boost::shared_ptr<const TypeBase> type_of_values_;
    unsigned int lower_bound_, upper_bound_;

public:
    template <class ValueType>
    Array(const ValueType &type, unsigned int min_size=0, unsigned int max_size=std::numeric_limits<unsigned int>::max() )
    : lower_bound_(min_size), upper_bound_(max_size)
    {
        type_of_values_ = boost::make_shared<ValueType>(type);
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

};

class Scalar : public TypeBase {
public:

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const {
        if (extensive) return stream;
        stream << "String (generic)";
        return stream;
    }
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
class Selection : public Scalar {
public:
    Selection();
    Selection(const std::vector<pair<Enum, std::string> > & enum_list);

    void add_value(const Enum value, const std::string &name);

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const {
        if (extensive) return stream;
        return stream;
    }

private:

};

class Bool : public Scalar {
public:
    Bool()
    {}

    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0)  const {
        stream << "Bool";
    }
};

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
 * Just for consistency, but is esentialy same as Scalar.
 */
class String : public Scalar {
};

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
