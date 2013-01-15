/*
 * accessors.hh
 *
 *  Created on: Mar 29, 2012
 *      Author: jb
 *
 *
 *

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
#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/static_assert.hpp>

#include "system/system.hh"
#include "system/exceptions.hh"

#include "input/input_type.hh"
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


/**
 * Class that represents base type of all enum types. We need it to return integer from a Selection input withou
 * knowing exact enum type. This class contains int and is convertible to int.
 */
class Enum {
public:
    Enum() : val_(0) {}
    Enum(int v) :val_(v) {}
    operator int() {return val_;}
    operator unsigned int() {return val_;}
private:
    int val_;
};


// Forward declaration
class IteratorBase;
template <class T> class Iterator;

/**
 * Class for storing and formating input address of an accessor (necessary for input errors detected after readed).
 *
 * TODO:
 * - store root TypeBase + array of iterator indices into Storage + cast from Abstract to Record
 * - at every level needs to copy internal vector, can be prohibitive only for
 *   large array of records ( one may need FastRecord accessor )
 * - possibly one can explicitly drop passing the address in constructor of an accessor/iterator
 *
 */
class InputAddress {
public:
private:
};

/**
 * @brief Accessor to the data with type \p Type::Record.
 *
 * This class provides access to the data through names -- key of the data fields.
 * It merge information from a @p Type::Record object, which describes valid keys and their types,
 * and reference into a StorageBase that provides access to actual values.
 *
 * The keys that are obligatory or has specified default value can be read either by the method template \p val<OutputType>
 * which returns the value of type \p OutputType that should be compatible with declared type (i.e. you can get unsigned int from
 * an input value of type Type:Integer, but not from an input value of type Type::String).
 *
 * The keys which are optional and has no default value you has to use
 * the method template \p find<OutputType> that returns an iterator to the type OutputType, which is invalid if the value
 * is missing on the input and no the default string is provided.
 *
 * Usage:
 @code
     using namespace Input;
     Record record = some_other_record.val<Record>("output_format");
     // reading an obligatory key or key with default value
     int n_digis = record.val<int>("number_of_substances");
     // reading an optional key
     Iterator<Array> it = record.find<Array>("substances_names");
     if ( it ) {
         it->copy_to( list_of_names );
     } else {
         // generic names of substances
     }

 @endcode
 *
 * @ingroup input_accessors
 *
 */
class Record {

public:
    /**
     * Default constructor creates an accessor to an empty storage.
     */
    Record()
    : record_type_(), storage_( NULL )
    {
    }



    /**
     * Copy constructor.
     */
    Record(const Record &rec)
    : record_type_(rec.record_type_), storage_(rec.storage_)
    {}


    /**
     * Constructs the accessor providing pointer \p store to storage node with list of data of the record and
     * type specification of the record given by parameter \p type.
     */
    Record(const StorageBase *store, const Type::Record type);

    /**
     * Returns value of given @p key if the declared key type (descendant of @p Input:Type:TypeBase) is convertible to the C++
     * class type given as the template parameter. If the key has no defined value
     * (either from input or a declared default value) it throws an exception. It throws also, if the
     * declared type do not match desired C++ type.
     *
     * This method can be used only for keys which are obligatory or has default value given at declaration.
     * The optional keys must use method @p find. Keys with default value at read time must use
     * the overloaded variant of method @p val or the method @p find
     *
     */
    template <class Ret>
    inline const Ret val(const string &key) const;

    /**
     * Same as the previous, but you can specify default value @p default_val that is used if the key is not specified at the input.
     * This method can be used only for keys declared with Default::reat_time().
     */
    template <class Ret>
    inline const Ret val(const string &key, const Ret default_val) const;

    /**
     * Returns iterator to the key if it exists or NULL Iterator if it doesn't.
     * This method must be used for keys which are optional or has default value provided at read time.
     */
    template <class Ret>
    inline Iterator<Ret> find(const string &key) const;

    /**
     * Returns true if the accessor is empty (after default constructor).
     */
    inline bool is_empty() const
    { return (storage_ == NULL); }


private:
    /// Corresponding Type::Record object.
    Input::Type::Record record_type_ ;

    /// Pointer to the corresponding array storage object.
    const StorageBase *storage_;
};



/**
 * @brief Accessor to the polymorphic input data of a type given by an AbstracRecord object.
 *
 * Provides conversion operator to the Record accessor in ordred to behave in the same way, but
 * further it provides method \p type() that can be used to call constructor of the class corresponding to the
 * input data.
 *
 * @ingroup input_accessors
 */

class AbstractRecord {
public:
    /**
     * Default constructor creates an accessor to an empty storage.
     */
    AbstractRecord()
    : record_type_(), storage_( NULL )
    {}

    /**
     * Copy constructor.
     */
    AbstractRecord(const AbstractRecord &rec)
    : record_type_(rec.record_type_), storage_(rec.storage_)
    {}


    /**
     * Constructs the accessor providing pointer \p store to storage node with list of data of the record and
     * type specification of the record given by parameter \p type.
     */
    AbstractRecord(const StorageBase *store, const Type::AbstractRecord type);

    /**
     * Implicit conversion to the \p Input::Record accessor. You can use \p Input::AbstractRecord in the same
     * way as the \p Input::Record.
     */
    operator Record();


    /**
     * Returns particular type selected from input. You can use it to construct particular type.
     *
     * @code
     * if (abstract_record.type() == MyClass.get_input_type())
     *      my_class = new MyClass(abstract_record);        // here the conversion to Input::Record is used
     * @endcode
     */
    Input::Type::Record type();


private:
    /// Corresponding Type::AbstractRecord object.
    Input::Type::AbstractRecord record_type_ ;

    /// Pointer to the corresponding array storage object.
    const StorageBase *storage_;
};



/**
 * @brief Accessor to input data conforming to declared Array.
 *
 * There are two possible ways how to retrieve data from Array accessor. First, you can use generic
 * @p copy_to function to copy the data into a given container. Second, you can get an Iterator<Type>
 * and iterate through the Array.
 *
 * In either case correspondence between resulting type (i.e. type of elements of the container or type of the Iterator)
 * and the type of the data in the Array is checked only once.
 *
 * Example of usage:
 * @code
 * Input::Array decay_array = in_rec.val<Input::Array>("decays");       // get accessor to an array stored under key 'decays'
 * int size = decay_array.size();                                       // get size of the actual arrya in the input (possibly for allocation)
 *
 * for(Input::Iterator<Input::Record> it = decay_array.begin<Input::Record>()   // pass through the array, that is array of records
 *          ; it != dacay_array.end(); ++it) {
 *
 *     Input::Iterator<double> it_hl = it->find<double>("half_life");   // check existence of an optional key
 *     if (it_hl) {
 *          double hl = *it_hl;
 *     } else {
 *          // use some other value, or turn-off the decay
 *     }
 *     Input::Array products = it->val<Input::Array>("products");       // read an obligatory key, theat conatins an array
 *     // ... process the array 'products'
 * }
 * @endcode
 *
 * @ingroup input_accessors
 */
class Array {
public:
    /**
     * Default constructor, empty accessor.
     */
    Array()
    : array_type_(Type::Bool()), storage_( NULL )
    {}

    /**
     * Copy constructor.
     */
    Array(const Array &ar)
    : array_type_(ar.array_type_), storage_(ar.storage_)
    {}

    /**
     * Constructs the accessor providing pointer \p store to storage node with list of data of the record and
     * type specification of the record given by parameter \p type.
     */
    Array(const StorageBase *store, const Type::Array type);

   /**
    * Returns iterator to the first element of input array. The template parameter is C++ type you want to
    * read from the array. Only types supported by Input::Interface::Iterator can be used.
    */
   template <class ValueType>
   Iterator<ValueType> begin() const;

   /**
    * Returns end iterator common to all iterators inner types.
    */
   inline IteratorBase end() const;

   /**
    * Actual size of the input array.
    */
   inline unsigned int size() const;

   /**
    * Method to fill the given container @p out with data in the input Array.
    * The container has to have methods @p clear and @p push_back. The C++ type of the
    * values in the container has to be supported by Iterator<T>.
    */
   template <class Container>
   void copy_to(Container &out) const;

private:
    /// Corresponding Type::Array.
    Input::Type::Array array_type_ ;

    /// Pointer to the corresponding array storage object.
    const StorageBase *storage_;
};


/*TODO:
 * Fast variant of RecordRead for reading array of records of same type. Has sense only for buffered input
 * storage and large data.
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


namespace internal {

/**
 *  This is primary type dispatch template. For given type it defines type that will be read from the storage.
 *  i.e. short int, char, and int are all translated to int.
 *
 *  TODO: us boost type traits to have this dispatch complete.
 *  following could work, but then the second dispatch do not work and probably has to be
 *  also implemented by mpl. The problem is, that TD<short int>::OT is not 'int', but unprocessed
 *  mpl construct. I don't know how to force compiler to process it before using it for the second dispatch.
 *
 @code
    template<class T>
    struct TD {
        typedef typename
                boost::mpl::if_< boost::is_integral<T>, int,
                    boost::mpl::if_< boost::is_floating_point<T>, double,
                        T
                    >
                >::type OT;
    };
  @endcode
  */
    template<class T>
    struct TD {
        typedef T OT;
    };
    /**
     * Secondary type dispatch. For every intermediate C++ type that can be read from input we have to define
     * read function from a given storage and Input type i.e. descendant of Input::Type::TypeBase.
     */

    template<class T>
    struct TypeDispatch;
} // close namespace internal


/**
 * Base class of input Iterator<Type> template. Main reason is possibility to construct
 * invalid iterator without template parameter ( used in \p Array::end() )
 */
class IteratorBase {
public:

    /**
     * Constructor of iterator without type and dereference methods.
     */
    IteratorBase(const StorageBase *storage, const unsigned int index)
    : storage_(storage), index_(index)
    {}

    /// Comparison of two Iterators. Do no compare types only position in the storage
    inline bool operator == (const IteratorBase &that) const;

    inline bool operator != (const IteratorBase &that) const;

    /**
     * Implicit conversion to bool. Returns true if iterator points to non-null storage.
     */
    inline operator bool() const;

protected:
    const StorageBase *storage_;
    unsigned int index_;
};



/**
 * This class behaves like iterator to type @p T (the template parameter), but in fact it is
 * iterator into input storage and also into tree of declarations through Input::Type classes.
 *
 * This class provides only limited functionality of iterators, namely prefix advance operator ++(),
 * dereference operator * (), dereference operator -> (), comparison operators == and != and implicit conversion to
 * bool (false in the case of the 'null' iterator).
 *
 *
 * @ingroup input_accessors
 */
template <class T>
class Iterator : public IteratorBase {
public:
    /// Converts C++ type @p T (template parameter) to 'DispatchType' from smaller set of types.
    typedef typename internal::TD<T>::OT DispatchType;
    /**
     * For small set of C++ types and accessor classes Record, AbstractRecord, and Array,
     * returns type of value given by dereference of the iterator (just add const to C++ types).
     */
    typedef typename internal::TypeDispatch<DispatchType>::ReadType OutputType;
    /**
     * A descendant of Input::Type::TypeBase that is appropriate to template parameter @p T.
     */
    typedef typename internal::TypeDispatch<DispatchType>::InputType InputType;


    /**
     * Constructor with Type of data
     */
    Iterator(const Input::Type::TypeBase &type,const StorageBase *storage, const unsigned int index)
    : IteratorBase(storage, index), type_( type_check_and_convert(type))
    {}

    /// Prefix. Advance operator.
    inline Iterator<T> &operator ++ ();

    /**
     *  Dereference operator * ; Shouldn't we return type T, i.e. try to cast from OutputType to T ??
     */
    inline OutputType operator *() const;

    /**
     *  Dereference operator can be used only for iterators to accessors Record, AbstractRecord, and Array.
     *
     *  TODO: make this without memory allocation.
     */
    inline OutputType *operator ->() const;


private:
    /// Iterator is not default constructible.
    Iterator();

    /**
     * Check that Type::TypeBase reference is in fact object of InputType
     * and returns converted copy (note that Type declaration objects are only handles.
     */
    static InputType type_check_and_convert(const Input::Type::TypeBase &type);

    /// Input type declaration.
    InputType type_;

    /**
     * temporary for -> operator, necessary only for T == DispatchType == OutputType
     * == any accessor: Record, AbstractRecord, Array
     *
     * for other types this is just an int.
     */
    mutable typename internal::TypeDispatch<DispatchType>::TmpType   temporary_value_;
};





namespace internal {

/**
 *  Template specializations for primary type dispatch.
 */
template<> struct TD<short int> { typedef int OT; };
template<> struct TD<unsigned short int> { typedef int OT; };
template<> struct TD<unsigned int> { typedef int OT; };
template<> struct TD<char> { typedef int OT; };
template<> struct TD<float> { typedef double OT; };

/**
 *  Template specializations for secondary type dispatch.
 */

// generic implementation accepts only enum types
template< class T>
struct TypeDispatch {
    BOOST_STATIC_ASSERT( ( boost::is_enum<T>::value || boost::is_same<T, Enum>::value ) );
    //BOOST_STATIC_ASSERT_MSG( boost::is_enum<T>::value , "TypeDispatch not specialized for given type." );

    typedef T TmpType;

    typedef Input::Type::Selection InputType;
    typedef const TmpType ReadType;
    static inline ReadType value(const StorageBase *s, const InputType&) { return ReadType( s->get_int() ); }
};


template<>
struct TypeDispatch<int> {
    typedef Input::Type::Integer InputType;
    typedef const int ReadType;
    typedef int TmpType;
    static inline ReadType value(const StorageBase *s, const InputType&) { return s->get_int(); }
};

template<>
struct TypeDispatch<bool> {
    typedef Input::Type::Bool InputType;
    typedef const bool ReadType;
    typedef int TmpType;
    static inline ReadType value(const StorageBase *s, const InputType&) { return s->get_bool(); }
};

template<>
struct TypeDispatch<double> {
    typedef Input::Type::Double InputType;
    typedef const double ReadType;
    typedef int TmpType;
    static inline ReadType value(const StorageBase *s, const InputType&) { return s->get_double(); }
};


template<>
struct TypeDispatch<string> {
    typedef Input::Type::String InputType;
    typedef const string ReadType;
    typedef int TmpType;
    static inline ReadType value(const StorageBase *s, const InputType&) { return s->get_string(); }
};


template<>
struct TypeDispatch<AbstractRecord> {
    typedef Input::Type::AbstractRecord InputType;
    typedef AbstractRecord ReadType;
    typedef AbstractRecord TmpType;
    static inline ReadType value(const StorageBase *s, const InputType& t) { return AbstractRecord(s, t); }

};


template<>
struct TypeDispatch<Record> {
    typedef Input::Type::Record InputType;
    typedef Record ReadType;
    typedef Record TmpType;
    static inline ReadType value(const StorageBase *s, const InputType& t) { return Record(s, t); }

};


template<>
struct TypeDispatch<Array> {
    typedef Input::Type::Array InputType;
    typedef Array ReadType;
    typedef Array TmpType;
    static inline ReadType value(const StorageBase *s, const InputType& t) { return Array(s,t); }

};

template<>
struct TypeDispatch<FilePath> {
    typedef Input::Type::FileName InputType;
    typedef FilePath ReadType;
    typedef int TmpType;
    static inline ReadType value(const StorageBase *s, const InputType& t) { return FilePath(s->get_string(), t.get_file_type() ); }

};




} // closing namespace internal





} // closing namespace Input


// include implementation of templates and inline methods
#include "accessors_impl.hh"



#endif
