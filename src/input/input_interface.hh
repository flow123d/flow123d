/*
 * input_interface.hh
 *
 *  Created on: Mar 29, 2012
 *      Author: jb
 *
 *
 *  todo:
 *  - accessor for AbstractRecord
 *
 *  - rename: access.hh, RecordAccessor , ArrayAccessor, (otazka co s iteraotorem, je divne, kdyz
 *    visi v namespacu Input
 *  - in iterator make dereference method with given pointer to default value string
 *  - key<Type>(key, default_value), otazka, zda defaultni hodnoty z deklarace neprirazovat az tady
 *  - presun vhodnych casti do *.cc
 *  - dokumentace
 *  - implementace operator -> without alocation (shared_ptr), i.e. put Accesors into Iterators
 *    Create corresponding acessor at construction of the iterator.
 *
 *  !!! how to pass instance of descendant of TypeBase through EI -
 *  - can not pass it directly since TypeBase is not copyconstructable
 *  - can not use shared_ptr for same reason
 *  - can not use C pointers since the refered object can be temporary
 *  solutions:
 *   - consistently move TypeBase to Pimpl design
 *   - provide virtual function make_copy, that returns valid shared_ptr
 */

#ifndef INPUT_INTERFACE_HH_
#define INPUT_INTERFACE_HH_

#include <vector>
#include <string>
#include "system.hh"
#include "system/exceptions.hh"

#include <boost/type_traits.hpp>

#include "input_type.hh"
#include "input/storage.hh"


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


namespace Input {


// exceptions and error_info types

TYPEDEF_ERR_INFO( EI_InputType, boost::shared_ptr<Input::Type::TypeBase>);
TYPEDEF_ERR_INFO( EI_RequiredType, const string );
TYPEDEF_ERR_INFO( EI_CPPRequiredType, const string );
TYPEDEF_ERR_INFO( EI_KeyName, const string);
DECLARE_EXCEPTION( ExcTypeMismatch, << "In Input::Interface:\n"
        << " can not make iterator with type " << EI_RequiredType::qval << "to get C++ type: " << EI_CPPRequiredType::qval
        << "\n since the key " << EI_KeyName::qval << " or array was declared with Input::Type : ";
//        << *(EI_InputType::ref(_exc))
        );




using std::string;

typedef StorageBase Storage;


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
    Record(const Record &rec)
    : record_type_(rec.record_type_), storage_(rec.storage_)
    {}


    /**
     * The only public constructor. It reads the file @p file_name, check the data structure against structure of @p record_type and
     * creates internal Storage of the data that can be read through resulting object.
     *
     * It takes Storage as a parameter instead of input stream until we know how to construct the Storage using data_tree class.
     *
     */
    Record(const Storage *store, const Type::Record type)
    : record_type_(type), storage_(store)
    {}

    /**
     * Returns value of given @p key if the declared key type (descendant of @p Input:Type:TypeBase) is convertible to the C++
     * class type given as the template parameter. If the key has no defined value
     * (either from input or a declared default value) it throws an exception. It throws also, if the
     * declared type do not match desired C++ type.
     */
    template <class Ret>
    inline const Ret key(const string &key) const {
        try {
            Type::Record::KeyIter key_it = record_type_.key_iterator(key);
            return *(Iterator<Ret>( *(key_it->type_), storage_, key_it->key_index));
        }
        catch (ExcTypeMismatch & e) {
            e << EI_CPPRequiredType(typeid(Ret).name()) << EI_KeyName(key);
            throw e;
        }

    }

    /**
     * Same as previous, but the value is returned through reference. This can be less verbose since you
     * need not to specify the template parameter. Further this method do not throw if the key is
     * missing instead it returns false in such a case. On the other hand you can not chain these methods to
     * get deeper into the input tree.
     */
    template <class Ret>
    inline bool has_key(const string &key, Ret &save_var) const {
        try {
            Type::Record::KeyIter key_it;
            if ( record_type_.has_key_iterator(key, key_it) ) {
                 save_var = *(Iterator<Ret>( *(key_it->type_), storage_, key_it->key_index));
                 return true;
            } else {
                return false;
            }

        } catch (ExcTypeMismatch & e) {
            e << EI_CPPRequiredType(typeid(Ret).name()) << EI_KeyName(key);
            throw e;
        }
    }

    /**
     * Same as previous, but with default value given now at read time instead at time of declaration. You must use this variant if the key is declared with
     * @p Default type @p read_time.
     */
    /*
    template <class Ret>
    Ret &key(const string &key, Ret  &default_value) const {
        int key_index = record_type_.key_index(key);
        ASSERT( key_index);
    }*/

private:
    Input::Type::Record record_type_ ;
    const Storage *storage_;


};





class AbstractRecord {
public:
    AbstractRecord(const AbstractRecord &rec)
    : record_type_(rec.record_type_), storage_(rec.storage_)
    {}


    /**
     * The only public constructor. It reads the file @p file_name, check the data structure against structure of @p record_type and
     * creates internal Storage of the data that can be read through resulting object.
     *
     * It takes Storage as a parameter instead of input stream until we know how to construct the Storage using data_tree class.
     *
     */
    AbstractRecord(const Storage *store, const Type::AbstractRecord type)
    : record_type_(type), storage_(store)
    {}


    operator Record()
    { return Record(storage_,type()); }


    /**
     * Returns particular type selected from input. You can use it to construct particular type.
     *
     * @code
     * if (abstract_record.type() == MyClass.get_input_record()) my_class = new MyClass(abstract_record);
     * @endcode
     */
    Input::Type::Record type()
    {
        unsigned int type_id = storage_->get_item(0)->get_int();
        return record_type_.get_descendant(type_id);
    }


private:
    Input::Type::AbstractRecord record_type_ ;
    const Storage *storage_;

};


/**
 * Is meant to be base class for storage iterators templated by type they points to.
 */

class IteratorBase {
public:
    IteratorBase(const Storage *storage, const unsigned int index)
    : storage_(storage), index_(index)
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

    Array(const Array &ar)
    : array_type_(ar.array_type_), storage_(ar.storage_)
    {}

    Array(const Storage *store, const Type::Array type)
    : array_type_(type), storage_(store)
    {}

   /**
    * Returns iterator to the first element of input array. The template parameter is C++ type you want to
    * read from the array. Only types supported by Input::Interface::Iterator can be used.
    */
   template <class ValueType>
   Iterator<ValueType> begin() const {
       try {
           return Iterator<ValueType>(array_type_.get_sub_type(), storage_, 0);
       }
       catch (ExcTypeMismatch & e) {
           e << EI_CPPRequiredType(typeid(ValueType).name()) << EI_KeyName("begin()");
           throw e;
       }
   }

   /**
    * Returns end iterator common to all iterators inner types.
    */
   inline IteratorBase end() const {
       return IteratorBase(storage_,storage_->get_array_size());
   }

   /**
    * Actual size of the input array.
    */
   inline unsigned int size() const {
       return storage_->get_array_size();
   }


   /**
    * Method to fill a container @p out with data in the input Array.
    * The container has to have methods @p clear and @p push_back. The C++ type of the
    * values in the container has to be supported by Iterator<T>.
    */
   template <class Container>
   void copy_to(Container &out) const {
       out.clear();
       Iterator<typename Container::value_type> it = begin<typename Container::value_type>();

       for(;it != end(); ++ it) {
           out.push_back(*it);
       }
   }
private:
    /// Corresponding Type::Array.
    Input::Type::Array array_type_ ;

    /// Pointer to the corresponding array storage object.
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

class FastRecord {

};

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
    static inline ReadType value(const Storage *s, const InputType&) { return s->get_int(); }
};

template<>
struct TypeDispatch<bool> {
    typedef Input::Type::Bool InputType;
    typedef const bool ReadType;
    static inline ReadType value(const Storage *s, const InputType&) { return s->get_bool(); }
};

template<>
struct TypeDispatch<double> {
    typedef Input::Type::Double InputType;
    typedef const double ReadType;
    static inline ReadType value(const Storage *s, const InputType&) { return s->get_double(); }
};


template<>
struct TypeDispatch<string> {
    typedef Input::Type::String InputType;
    typedef const string ReadType;
    static inline ReadType value(const Storage *s, const InputType&) { return s->get_string(); }
};


template<>
struct TypeDispatch<AbstractRecord> {
    typedef Input::Type::AbstractRecord InputType;
    typedef AbstractRecord ReadType;
    static inline ReadType value(const Storage *s, const InputType& t) { return AbstractRecord(s, t); }

    static Record temporary_;
};


template<>
struct TypeDispatch<Record> {
    typedef Input::Type::Record InputType;
    typedef Record ReadType;
    static inline ReadType value(const Storage *s, const InputType& t) { return Record(s, t); }

    static Record temporary_;
};

template<>
struct TypeDispatch<Array> {
    typedef Input::Type::Array InputType;
    typedef Array ReadType;
    static inline ReadType value(const Storage *s, const InputType& t) { return Array(s,t); }

    static Array temporary_;
};

template<>
struct TypeDispatch<FilePath> {
    typedef Input::Type::FileName InputType;
    typedef FilePath ReadType;
    static inline ReadType value(const Storage *s, const InputType& t) { return FilePath(s->get_string(), t.get_file_type() ); }

    static Array temporary_;
};


/**
 * This class behaves like iterator to type @p T (the template parameter), but in fact it is
 * iterator into input storage and also into tree of declarations through Input::Type classes.
 *
 */
template <class T>
class Iterator : public IteratorBase {
public:
    // Appropriate Dispatch helper class that contains things that depends on T.
    typedef typename TD<T>::OT DispatchType;
    // Appropriate type actually given be dispatcher.
    typedef typename TypeDispatch<DispatchType>::ReadType OutputType;
    // Appropriate declaration type - descendant of Type::TypeBase
    typedef typename TypeDispatch<DispatchType>::InputType InputType;

    static InputType type_check_and_convert(const Input::Type::TypeBase &type)
    {
        if ( typeid(type) == typeid(InputType)) {
            return static_cast< const InputType & >( type );
        } else {
            THROW( ExcTypeMismatch()
                    //<< EI_InputType( boost::make_shared<Input::Type::TypeBase>(type) )
                    << EI_RequiredType( typeid(InputType).name() ) );
        }
    }

    /**
     * Constructor needs Type of data
     */
    Iterator(const Input::Type::TypeBase &type,const Storage *storage, const unsigned int index)
    : IteratorBase(storage, index), type_( type_check_and_convert(type))
    {}




    /// Prefix. Advance operator.
    inline Iterator<T> &operator ++ ()
    {
        index_++;
        return *this;
    }

    /**
     *  Dereference operator * ; shouldn;t we return type T, i.e. try to cast from OutputType to T ??
     */
    inline OutputType operator *() const
    {
        const Storage *s = storage_->get_item(index_);
        if (s && ! s->is_null()) {
            return TypeDispatch<DispatchType>::value(s, type_);
        } else {
            // error no value
        }
    }

    /**
     *  Dereference operator can be used only for iterators to Record or Array.
     */
    inline OutputType *operator ->() const
    {
        BOOST_STATIC_ASSERT( ( boost::is_base_of<Record, OutputType>::value ||
                               boost::is_base_of<Array, OutputType>::value  ) );

        // we have to make save temporary
        OutputType xx = this->operator*();
        return boost::make_shared<OutputType>( xx  ).get();

    }



private:
    /// Iterator is not default constructable.
    Iterator();


    InputType type_;


};

/*********************************************************************
 * Implementation of Array accessor
 */



/*********************************************************************
 * Implementation of Iterator<T>
 */

} // closing namespace Input

#endif
