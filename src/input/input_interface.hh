/*
 * input_interface.hh
 *
 *  Created on: Mar 29, 2012
 *      Author: jb
 *
 *
 *  todo:
 *  - in iterator make dereference method with given pointer to default value string
 *  - key<Type>(key, default_value), otazka, zda defaultni hodnoty z deklarace neprirazovat az tady
 *  - pouhy dotaz na existenci klice
 *  - presun vhodnych casti do *.cc
 *  - dokumentace
 *  - deklarovat tridu Path a mit key<Path>, ktera bude produkovat uplne cesty, varianta je
 *
 *
 */
#include <vector>
#include <string>
#include "system.hh"
#include "system/exceptions.hh"

#include <boost/type_traits.hpp>

#include "input_type.hh"

//#include "Generic_node.hpp"
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

//* proper Generic_node stub (all in stub is implemented with same data types)

class Generic_node {
public:
    const int get_int() const
        {return 0;}
    const double get_double() const
        {return 0.0;}
    const bool get_bool() const
        {return false;}
    const string get_string() const
        {return *( new string(""));} // memory leak

    const Generic_node &get_item(const size_t index) const
        {return *( new Generic_node() );} // memory leak
    bool not_null() const { return true;}
    bool is_null() const { return false;}
    size_t get_array_size() const {return 1;}
};

namespace Input {
namespace Interface {

// exceptions and error_info types

struct TypeMismatchExcp : virtual FlowException {};
TYPEDEF_ERR_INFO( InputType, const Input::Type::TypeBase *);



using std::string;

typedef Generic_node Storage;

enum ErrorCode {
    no_error = 0,       ///< no error occured
    no_value = 1        ///< no value on input and no default value provided
};



template <class T>
class Iterator;



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
     * It takes Storage as a parameter instead of input stream until we know how to construct the Storage using data_tree class.
     *
     */
    Record(const Storage &store, const Type::Record &type)
    : record_type_(&type), storage_(&store)
    {}

    /**
     * Returns value of given @p key if the declared key type (descendant of @p Input:Type:TypeBase) is convertible to the C++
     * class type given as the template parameter. If the key has no defined value (either from input or a declared default value) the error is reported.
     */
    template <class Ret>
    inline const Ret key(const string &key) const {
        try {
            Type::Record::KeyIter key_it = record_type_->key_iterator(key);
            return *(Iterator<Ret>( *(key_it->type_), *storage_, key_it->key_index));
        }
        catch (FlowException & e) {
            const Input::Type::TypeBase * key_type = *( boost::get_error_info<InputType_EI>(e) );
            std::stringstream s_key;
            s_key << *key_type;
            std::stringstream s_record;
            s_record << *record_type_;
            xprintf(Err, "Error: Can not return value of C++ type '%s' from key '%s' of type '%s' in record type:\n %s",
                    typeid(Ret).name(), key.c_str(), s_key.str().c_str(), s_record.str().c_str());
        }

    }

    /**
     * Same as previous, but if the key has no value (or on other error) an undefined value is returned with nonzero error @p code.
     * This can be used to read optional keys.
     */
    /*
    template <class Ret>
    Ret &key(const string &key, ErrorCode &code) const {
        int key_index = record_type_->key_index(key);
        ASSERT( key_index);
    }

    /**
     * Same as previous, but with default value given now at read time instead at time of declaration. You must use this variant if the key is declared with
     * @p DefaultValue type @p read_time.
     */
    /*
    template <class Ret>
    Ret &key(const string &key, Ret  &default_value) const {
        int key_index = record_type_.key_index(key);
        ASSERT( key_index);
    }*/

private:
    const Input::Type::Record *record_type_ ;
    const Storage *storage_;


};

class IteratorBase;
/**
 * Is meant to be base class for storage iterators templated by type they points to.
 */

class IteratorBase {
public:
    IteratorBase(const Storage &storage, const unsigned int index)
    : storage_(&storage), index_(index)
    {}
    /// Comparison of two Iterators.
    inline bool operator == (const IteratorBase &that) const
            { return ( storage_  == that.storage_  && index_ == that.index_); }

    inline bool operator != (const IteratorBase &that) const
            { return ! ( *this == that ); }
protected:
    const Storage *storage_;
    unsigned int index_;
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
    Array(const Storage &store, const Type::Array &type)
    : array_type_(&type), storage_(&store)
    {}

   /**
    * Returns iterator to the first element of input array. The template parameter is C++ type you want to
    * read from the array. Only types supported by Input::Interface::Iterator can be used.
    */
   template <class ValueType>
   Iterator<ValueType> begin() {
       try {
           return Iterator<ValueType>(array_type_->get_sub_type(), *storage_, 0);
       }
       catch (TypeMismatchExcp & e) {
           const Input::Type::TypeBase * key_type = *( boost::get_error_info<InputType_EI>(e) );
           std::stringstream s_key;
           s_key << *key_type;

           xprintf(Err, "Error: Can not get iterator pointing to C++ type '%s' from array of values of type:\n %s\n",
                   typeid(ValueType).name(), s_key.str().c_str());
       }

   }

   /**
    * Returns end iterator common to all iterators inner types.
    */
   inline IteratorBase end() {
       return IteratorBase(*storage_,storage_->get_array_size());
   }

   /**
    * Method to fill a container @p out with data in the input Array. The container has to have methods @p clear and @p push_back. The C++ type of the
    * values in the container has to be supported by Iterator<T>.
    */
   template <class Container>
   void copy_to(Container &out) {
       out.clear();
       Iterator<typename Container::value_type> it = begin<typename Container::value_type>();

       for(;it != end(); ++ it) {
           out.push_back(*it);
       }
   }
private:
    const Input::Type::Array *array_type_ ;
    const Storage *storage_;
};


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



class Record;
class Array;

/**
 *  This is primary type dispatch template. OT has to be convertible to T.
 */
template<class T>
struct TD {
    typedef T OT;
};

template<>
struct TD<short int> {
    typedef int OT;
};

template<>
struct TD<char> {
    typedef int OT;
};

template<>
struct TD<float> {
    typedef double OT;
};

/**
 * Secondary type dispatch. For every intermediate C++ type that can be read from input we have to define
 * read function from a given storage and Input type i.e. descendant of Input::Type::TypeBase.
 */
template<class T>
struct TypeDispatch {

};

template<>
struct TypeDispatch<int> {
    typedef Input::Type::Integer InputType;
    typedef const int ReadType;
    static inline ReadType value(const Storage &s, const InputType&) { return s.get_int(); }
};

template<>
struct TypeDispatch<bool> {
    typedef Input::Type::Bool InputType;
    typedef const bool ReadType;
    static inline ReadType value(const Storage &s, const InputType&) { return s.get_bool(); }
};

template<>
struct TypeDispatch<double> {
    typedef Input::Type::Double InputType;
    typedef const double ReadType;
    static inline ReadType value(const Storage &s, const InputType&) { return s.get_double(); }
};


template<>
struct TypeDispatch<string> {
    typedef Input::Type::String InputType;
    typedef const string ReadType;
    static inline ReadType value(const Storage &s, const InputType&) { return s.get_string(); }
};

template<>
struct TypeDispatch<Record> {
    typedef Input::Type::Record InputType;
    typedef Record ReadType;
    static inline ReadType value(const Storage &s, const InputType& t) { return Record(s, t); }
};

template<>
struct TypeDispatch<Array> {
    typedef Input::Type::Array InputType;
    typedef Array ReadType;
    static inline ReadType value(const Storage &s, const InputType& t) { return Array(s,t); }
};




/**
 *
 */
template <class T>
class Iterator : public IteratorBase {
public:
    typedef typename TD<T>::OT DispatchType;
    typedef typename TypeDispatch<DispatchType>::ReadType OutputType;

    /**
     * Constructor needs Type of data
     */
    Iterator(const Input::Type::TypeBase &type,const Storage &storage, const unsigned int index)
    : IteratorBase(storage, index)
    {
        if (typeid(type) == typeid(InputType)) {
            type_ = static_cast< const InputType * >( &type );
        } else {
            //DBGMSG("type: '%s' input type: '%s'\n", typeid(type).name(), typeid(InputType).name() );
            throw TypeMismatchExcp() << InputType_EI(&type) ;
        }
    }



    /// Prefix. Advance operator.
    inline Iterator<T> &operator ++ ()
    {
        index_++;
        return *this;
    }

    ///  * dereference operator
    inline OutputType operator *() const
    {
        const Storage &s = storage_->get_item(index_);
        if (s.not_null()) {
            return TypeDispatch<DispatchType>::value(s, *(type_));
        } else {
            // error no value
        }
    }

    /*
     *  dereference operator should be implemented only for Record
    inline const T * operator ->() const
    {
        return &(*(*this));
    }
    */



private:
    /// Iterator is not default constructable.
    Iterator();

    typedef typename TypeDispatch<DispatchType>::InputType InputType;
    const InputType *type_;


};




} // closing nemaspace Read
} // closing namespace Input

