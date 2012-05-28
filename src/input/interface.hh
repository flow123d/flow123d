/*
 * interface.hh
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



namespace Input {

using std::string;

// exceptions and error_info types

// throwed in Iterator<>
TYPEDEF_ERR_INFO( EI_InputType, const string);
TYPEDEF_ERR_INFO( EI_RequiredType, const string );
TYPEDEF_ERR_INFO( EI_CPPRequiredType, const string );
TYPEDEF_ERR_INFO( EI_KeyName, const string);
DECLARE_EXCEPTION( ExcTypeMismatch, << "Key:" << EI_KeyName::qval
        << ". Can not construct Iterator<T> with C++ type T=" << EI_CPPRequiredType::qval << ";\n"
        << "can not convert Type: " << EI_InputType::qval << " to: " << EI_RequiredType::qval
        );

// throwed in Record, Array, AbstractRecord
TYPEDEF_ERR_INFO( EI_AccessorName, const string );
DECLARE_EXCEPTION( ExcAccessorForNullStorage, << "Can not create " << EI_AccessorName::val << " from StorageNull.");




template <class T>
class Iterator;


/**
 * @brief Accessor to the data conforming to a declared record.
 *
 * This class provides access to the data through string names -- key of the data fields.
 * It merge information form a @p Type::Record object, which describes valid keys and their types,
 * and reference into a StorageBase that provides access to actual values.
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
     * creates internal StorageBase of the data that can be read through resulting object.
     *
     * It takes StorageBase as a parameter instead of input stream until we know how to construct the StorageBase using data_tree class.
     *
     */
    Record(const StorageBase *store, const Type::Record type)
    : record_type_(type), storage_(store)
    {
        if (store->is_null())
            THROW( ExcAccessorForNullStorage() << EI_AccessorName("Record") );
    }

    /**
     * Returns value of given @p key if the declared key type (descendant of @p Input:Type:TypeBase) is convertible to the C++
     * class type given as the template parameter. If the key has no defined value
     * (either from input or a declared default value) it throws an exception. It throws also, if the
     * declared type do not match desired C++ type.
     */
    template <class Ret>
    inline const Ret val(const string &key) const {
        try {
            Type::Record::KeyIter key_it = record_type_.key_iterator(key);

            ASSERT(! key_it->default_.is_optional(),
                    "The key %s is declared as optional, you have to use Record::find instead.\n", key.c_str());

            Iterator<Ret> it = Iterator<Ret>( *(key_it->type_), storage_, key_it->key_index);
            return *it;
        }
        // we catch all possible exceptions
        catch (Type::Record::ExcRecordKeyNotFound & e) {
            throw;
        }
        catch (ExcTypeMismatch & e) {
            e << EI_CPPRequiredType(typeid(Ret).name()) << EI_KeyName(key);
            throw;
        }
        catch (ExcStorageTypeMismatch &e) {
            throw;
        }
        catch (ExcAccessorForNullStorage &e) {
            throw;
        }
    }



    template <class Ret>
    inline Iterator<Ret> find(const string &key) const {
        try {
            Type::Record::KeyIter key_it = record_type_.key_iterator(key);
            return Iterator<Ret>( *(key_it->type_), storage_, key_it->key_index);
        }
        // we catch all possible exceptions
        catch (Type::Record::ExcRecordKeyNotFound & e) {
            throw;
        }
        catch (ExcTypeMismatch & e) {
            e << EI_CPPRequiredType(typeid(Ret).name()) << EI_KeyName(key);
            throw;
        }

    }


private:
    Input::Type::Record record_type_ ;
    const StorageBase *storage_;


};





class AbstractRecord {
public:
    AbstractRecord(const AbstractRecord &rec)
    : record_type_(rec.record_type_), storage_(rec.storage_)
    {}


    /**
     * The only public constructor. It reads the file @p file_name, check the data structure against structure of @p record_type and
     * creates internal StorageBase of the data that can be read through resulting object.
     *
     * It takes StorageBase as a parameter instead of input stream until we know how to construct the StorageBase using data_tree class.
     *
     */
    AbstractRecord(const StorageBase *store, const Type::AbstractRecord type)
    : record_type_(type), storage_(store)
    {
        if (store->is_null())
            THROW( ExcAccessorForNullStorage() << EI_AccessorName("AbstractRecord") );
    }


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
    const StorageBase *storage_;

};


/**
 * Is meant to be base class for storage iterators templated by type they points to.
 */

class IteratorBase {
public:



    IteratorBase(const StorageBase *storage, const unsigned int index)
    : storage_(storage), index_(index)
    {}
    /// Comparison of two Iterators.
    inline bool operator == (const IteratorBase &that) const
            { return ( storage_  == that.storage_  && index_ == that.index_); }

    inline bool operator != (const IteratorBase &that) const
            { return ! ( *this == that ); }
protected:
    const StorageBase *storage_;
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

    Array(const StorageBase *store, const Type::Array type)
    : array_type_(type), storage_(store)
    {
        if (store->is_null())
            THROW( ExcAccessorForNullStorage() << EI_AccessorName("Array") );
    }

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
    const StorageBase *storage_;
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
    static inline ReadType value(const StorageBase *s, const InputType&) { return s->get_int(); }
};

template<>
struct TypeDispatch<bool> {
    typedef Input::Type::Bool InputType;
    typedef const bool ReadType;
    static inline ReadType value(const StorageBase *s, const InputType&) { return s->get_bool(); }
};

template<>
struct TypeDispatch<double> {
    typedef Input::Type::Double InputType;
    typedef const double ReadType;
    static inline ReadType value(const StorageBase *s, const InputType&) { return s->get_double(); }
};


template<>
struct TypeDispatch<string> {
    typedef Input::Type::String InputType;
    typedef const string ReadType;
    static inline ReadType value(const StorageBase *s, const InputType&) { return s->get_string(); }
};


template<>
struct TypeDispatch<AbstractRecord> {
    typedef Input::Type::AbstractRecord InputType;
    typedef AbstractRecord ReadType;
    static inline ReadType value(const StorageBase *s, const InputType& t) { return AbstractRecord(s, t); }

};


template<>
struct TypeDispatch<Record> {
    typedef Input::Type::Record InputType;
    typedef Record ReadType;
    static inline ReadType value(const StorageBase *s, const InputType& t) { return Record(s, t); }

};

template<>
struct TypeDispatch<Array> {
    typedef Input::Type::Array InputType;
    typedef Array ReadType;
    static inline ReadType value(const StorageBase *s, const InputType& t) { return Array(s,t); }

};

template<>
struct TypeDispatch<FilePath> {
    typedef Input::Type::FileName InputType;
    typedef FilePath ReadType;
    static inline ReadType value(const StorageBase *s, const InputType& t) { return FilePath(s->get_string(), t.get_file_type() ); }

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



    /**
     * Constructor needs Type of data
     */
    Iterator(const Input::Type::TypeBase &type,const StorageBase *storage, const unsigned int index)
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
        const StorageBase *s = storage_->get_item(index_);

        ASSERT( s, "NULL pointer in storage!!! \n");

        return TypeDispatch<DispatchType>::value(s, type_);
    }

    /**
     *  Dereference operator can be used only for iterators to Record or Array.
     */
    inline OutputType *operator ->() const
    {
        BOOST_STATIC_ASSERT( ( boost::is_same<Record, OutputType>::value ||
                               boost::is_same<AbstractRecord, OutputType>::value ||
                               boost::is_same<Array, OutputType>::value  ) );

        // we have to make save temporary
        OutputType xx = this->operator*();
        return boost::make_shared<OutputType>( xx  ).get();

    }

    /**
     * Implicit conversion to bool. Returns true if iterator points to non-null storage.
     */
    inline operator bool() {
        const StorageBase *s = storage_->get_item(index_);
        return ( s && ! s->is_null() );
    }


private:
    /// Iterator is not default constructible.
    Iterator();

    /**
     * Check that Type::TypeBase reference is in fact object of InputType
     * and returns converted copy (note that Type declaration objects are only handles.
     */
    static InputType type_check_and_convert(const Input::Type::TypeBase &type)
    {
        if ( typeid(type) == typeid(InputType)) {
            return static_cast< const InputType & >( type );
        } else {
            THROW( ExcTypeMismatch()
                    << EI_InputType( type.type_name() )
                    << EI_RequiredType( typeid(InputType).name() ) );
        }
    }

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
