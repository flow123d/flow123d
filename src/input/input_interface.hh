/*
 * input_interface.hh
 *
 *  Created on: Mar 29, 2012
 *      Author: jb
 */
#include <vector>
#include <string>
#include "system.hh"
#include "input_type.hh"

#include <boost/type_traits.hpp>

namespace Input {
namespace Interface {

using std::string;


enum ErrorCode {
    none = 0,
    key_not_found=1
};
class Storage;

/**
 * @brief Accessor to the data conforming to a declared record.
 *
 * This class provides access to the data through string names -- key of the data fields.
 * It merge information form a @p Type::Record object, which describes valid keys and their types,
 * and reference into a Storage that provides access to actual values.
 *
 *
 * Usage:
 * Read::Record record= ... from somewhere
 * int i = record.key<int>("integer_key");
 * record = record.key<Read::Record>("sub_record");
 *
 *
 */
class Record {

public:
    /**
     * The only public constructor. It reads the file @p file_name, check the data structure against structure of @p record_type and
     * creates internal Storage of the data that can be read through resulting object.
     *
     */
    Record(Storage * store, boost::shared_ptr< Type::Record > &type)
    : record_type_(type), storage_(store)
    {}

    /**
     * Returns value of given @p key if the declared key type (descendant of @p Input:Type:TypeBase) is convertible to the C++
     * class type given as the template parameter. If the key has no defined value (either from input or a declared default value) the error is reported.
     */
    template <class Ret>
    inline Ret &key(const string &key) const {
        int key_index = record_type_->key_index(key);
        ASSERT( key_index);
    }

    /**
     * Same as previous, but if the key has no value (or on other error) an undefined value is returned with nonzero error @p code.
     * This can be used to read optional keys.
     */
    template <class Ret>
    inline Ret &key(const string &key, ErrorCode &code) const {
        int key_index = record_type_->key_index(key);
        ASSERT( key_index);
    }

    /**
     * Same as previous, but with default value given now at read time instead at time of declaration. You must use this variant if the key is declared with
     * @p DefaultValue type @p read_time.
     */
    template <class Ret>
    inline Ret &key(const string &key, Ret  &default_value) const {
        int key_index = record_type_->key_index(key);
        ASSERT( key_index);
    }

private:
    boost::shared_ptr<const ::Input::Type::Record> record_type_;
    Storage * storage_;


};
/*
template <class T>
inline T get_key(TypeBase *t, Storage *s) {

    if (boost::is_integral)::value == true) return get_key_int<T>(t,s);
    else if (boost::is_enum<T>::value == true) return get_key_enum<T>(t, s);
    else if (boost::is_floating_point<T>::value == true) return get_key_double<T>(t,s);
    else if (boost::is_same<T, std::string>::value == true) return get_key_string<T>(t,s);
    else if (boost::is_same<T, Read::Record>::value == true) return get_key_record<T>(t,s);
    else if (boost::is_same<T, Read::Array>::value == true) return get_key_array<T>(t,s);
    else return get_key_error<T>(t,s);
}

template <class T>
inline T get_key_error(TypeBase *t, Storage *s) {
    // this function is defined only for some types
    BOOST_STATIC_ASSERT(true);
}

template <T>
T get_key_int<T>(TypeBase *t, Storage *s) {
    if (typeid(t) == typeid(Type::Integer)) {
        return s.get_int();
    }
}

/**
 * Fast variant of RecordRead for reading array of records of same type.
 *
 * usage:
 *
 * Array<FastRecordReader> a_of_fr;
 * FastRecordReader::iterator<double> iter_x = a_of_fr.get_type().iter_of_key<double>("x");
 * FastRecordReader::iterator<double> iter_y = a_of_fr.get_type().iterof_key<double>("y");
 * for(Array<FastRecordReader>::iterator it=a_of_fr.begin(); it != a_of_fr.end(); ++it) {
 *      it->fast_get(iter_x);
 *      it->fast_get(iter_y);
 * }
 *
 */





class IteratorBase {
public:
/*
protected:
    template<>
    Read::Array &get_value<Read::Array>() {
        return Array( storage_->item(index_)->get_storage() );
    }

    template<>
    Read::Record &get_value<Read::Array>() {
        return Record( storage_->item(index_)->get_storage() );
    }


    template<>
    T &get_value<string>() {
        return storage_->item(index_)->get_string();
    }

    template<>
    T &get_value<int>() {
        return storage_->item(index_)->get_int();
    }
    template<>
    T &get_value<double>() {
        return storage_->item(index_)->get_double();
    }
    template<>
    T &get_value<bool>() {
        return storage_->item(index_)->get_bool();
    }
*/
    Storage * storage_;
    int index_;
};

template <class T>
class Iterator : public IteratorBase {
public:
    ///  * dereference operator
    inline const T & operator *() const
    { //return get_value<T>();
    }

    inline const T * operator ->() const
    {
        return &(*(*this));
    }

    /// Comparison of two FullIterator.
    inline bool operator == (const IteratorBase &that) const
            { return ( this->storage_  == that.storage_  && this->index_ == that.index_); }

    inline bool operator != (const IteratorBase &that) const
            { return ! ( *this == that ); }

    /// Prefix. Advance operator.
    inline Iterator<T> &operator ++ ()
    {
        index_++;
        return *this;
    }

private:
};

/**
 * @brief Accessor to input data conforming to declared Array.
 *
 * There are two possible ways how to retrieve data from Array accessor. First, you can use generic
 * @p copy_to function to copy the data into a given container. Second, you can get an Iterator<Type>
 * and pass through the Array
 *
 * In either case correspondence between resulting type (i.e. type of elements of the container or type of the Iterator)
 * and the type of the data in the Array is checked only once.
 */

class Array {
public:
   /**
    *
    */
   template <class ValueType>
   Iterator<ValueType> begin() {

   }

   IteratorBase end() {

   }
   /**
    *
    */
   template <class Container>
   void copy_to(Container &out) {

   }
};

/**
 * Object for storing data tree. One node can hold:
 * - int, double, bool
 * - pointer to string
 * - pointer to array of storages
 * - special state: NULL (no data), REF (reference to other place of storage tree - maybe in future, until then we copy everything),
 *   INCLUDE (have to read another file to provide the value, this may be possible only through particular readers)
 *   ...
 *
 * Not all readers has to use Storage for accessing the input data !!
 */
class Storage {

};

} // closing nemaspace Read
} // closing namespace Input

