/*
 * InputType.hh
 *
 *  Created on: Mar 28, 2012
 *      Author: jb
 */

#ifndef INPUTTYPE_HH_
#define INPUTTYPE_HH_

#include <limits>
#include <ios>
#include <map>
#include <string>

#include "system.hh"

namespace Input {

#include <const_hashes.h>

using std::string;

enum DefaultType {
    no_default,
    default_at_read,
    default_at_declaration
};

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
 */
class TypeBase {
protected:

public:
    /**
     * Fundamental method for output documentation of Input Types.
     */
    virtual std::ostream& documentation(std::ostream& stream) const = 0;
};

/**
 * For convenience we provide redirection operator for output documentation of Input:Type classes.
 */
std::ostream& operator<<(std::ostream& stream, const TypeBase& type) {
    return type.documentation(stream);
}

class AbstractRecord : public TypeBase {

};

class Scalar;
class Array;

/**
 * Possible types of default values:
 * - given at declaration time
 * - provided at read time (not part of Record specification)
 * - none
 */
class Record : public TypeBase {
private:
    struct Key {
        const string &key_;
        const string &description_;
        const TypeBase &type_;
        enum DefaultType dflt_type_;
        const string &default_;
    };

    std::map<unsigned int, Key> keys;

    unsigned int hash(const string &str) {
        return ( CONSTHASH(str.c_str() ) );
    }
public:
    void declare_key(const string &key,const Scalar &type,const string &description,
                     const string &default_value )
    {
        declare_key_(key, (TypeBase &) type,description,default_at_declaration, default_value, hash(key) );
    }

    void declare_key(const string &key,const Scalar &type,const string &description,
                     const DefaultType dflt_type)

    {
        if (dflt_type != default_at_declaration)
            declare_key_(key,(TypeBase &)type,description,dflt_type, "", hash(key) );
        else
            xprintf(Err, "Can not call declare_key with 'default_at_declaration' without providing default value.\n");
    }

    void declare_key(const string &key,const TypeBase &type,const string &description)
    {
        declare_key_(key,type,description,no_default, "", hash(key) );
    }

    /**
     * The most general declaration of a key.
     */
    void declare_key_(const string& key, const TypeBase& type, const string& description,
                     const DefaultType dflt_type, const string& dflt, const unsigned int key_hash)
    {
        const Key tmp_key= {key, description, type, dflt_type, dflt};
        keys.insert( std::make_pair(key_hash, tmp_key) );
    }

};


class Array : public TypeBase {
private:
    const TypeBase &type_of_values_;
    unsigned int lower_bound_, upper_bound_;

public:
    Array(const TypeBase &type, unsigned int min_size=0, unsigned int max_size=std::numeric_limits<unsigned int>::max() )
    : type_of_values_(type),lower_bound_(min_size), upper_bound_(max_size)
    {}

    virtual void documentation(std::ostream &stream) {
        stream << "Array, limits: [" << lower_bound_ << ", " << upper_bound_ << "], values of type:\n";
        stream << "       ";
        type_of_values_.documentation(stream);
    }
};

class Scalar : public TypeBase {

};

class ScalarInteger : public Scalar {
private:
    int lower_bound_, upper_bound_;
public:
    ScalarInteger(int lower_bound=std::numeric_limits<int>::min(), int upper_bound=std::numeric_limits<int>::max())
    : lower_bound_(lower_bound), upper_bound_(upper_bound_)
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

    virtual void documentation(std::ostream &stream) {
        stream << "Integer in [" << lower_bound_ << ", " << upper_bound_ << "]";
    }

};

class ScalarDouble : public Scalar {
private:
    double lower_bound_, upper_bound_;
public:
    ScalarDouble(double lower_bound=std::numeric_limits<double>::min(), double upper_bound=std::numeric_limits<double>::max())
    : lower_bound_(lower_bound), upper_bound_(upper_bound_)
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

    virtual void documentation(std::ostream &stream) {
        stream << "Double in [" << lower_bound_ << ", " << upper_bound_ << "]";
    }
};


} // closing namespace InputType


#endif /* INPUTTYPE_HH_ */
