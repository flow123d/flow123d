/*
 * accessors.hh
 *
 *  Created on: Mar 29, 2012
 *      Author: jb
 *
 *
 *  TODO:
 *  - decide which part of interface has to be optimized ( probably nothing until we
 *    implement reader for HDF5, XML or large Raw data files, and try to use the same input interface for input of large data)
 *  - then make inlined only neccessary functions  and carefully move as much as possible into accessors.cc including explicit instantiation of
 *    support classes. This should speedup compilation of the code that use the accessors.
 *
 *  - implement operator -> without allocation (shared_ptr), i.e. put Accesors into Iterators
 *    Create corresponding accessor at construction of the iterator.
 *
 */

#ifndef INPUT_INTERFACE_HH_
#define INPUT_INTERFACE_HH_

#include <vector>
#include <string>
#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/static_assert.hpp>
#include <boost/shared_ptr.hpp>

#include "system/system.hh"
#include "system/exceptions.hh"

#include "input/input_type.hh"
#include "input/factory.hh"
#include "input/storage.hh"







namespace Input {

using std::string;


/**
 * @brief Base of exceptions due to user input.
 *
 * Base class for "input exceptions" that are exceptions caused by incorrect input from the user
 * not by an internal error.
 *
 * @ingroup exceptions
 */
class Exception : public virtual ExceptionBase
{
public:
    const char * what () const throw ();
    virtual ~Exception() throw () {};
};


/**
 *  Declaration of error info class for passing Input::Address through exceptions.
 *  Is returned by input accessors : Input::Record, Input::Array, etc.
 *
 *  TODO: if Address class is persistent (every copy is self contented, we can use Address instead of std::string.
 *  see also ei_address methods.
 */
TYPEDEF_ERR_INFO( EI_Address, const std::string);


/**
 * @brief Macro for simple definition of input exceptions.
 *
 * Works in the same way as @p DECLARE_EXCEPTION, just define class derived from
 * @p InputException. Meant to be used for exceptions due to wrong input from user.
 *
 * @ingroup exceptions
 */
#define DECLARE_INPUT_EXCEPTION( ExcName, Format)                             \
struct ExcName : public virtual ::Input::Exception {                          \
     virtual void print_info(std::ostringstream &out) const {                 \
         using namespace internal;                                            \
         ::internal::ExcStream estream(out, *this);                           \
         estream Format														  \
      	  	  	  << "\nAt input address: " 			   					\
      	  	  	  << ::Input::EI_Address::val; 								\
      	 out << std::endl;													\
     }                                                                      \
     virtual ~ExcName() throw () {}                                         \
}

/**
 * Simple input exception that accepts just string message.
 */
DECLARE_INPUT_EXCEPTION(ExcInputMessage, << EI_Message::val );





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

// throwed in Address
TYPEDEF_ERR_INFO( EI_ParamName, const string);
DECLARE_EXCEPTION( ExcAddressNullPointer, << "NULL pointer in " << EI_ParamName::val << " parameter.");




/**
 * Class that works as base type of all enum types. We need it to return integer from a Selection input without
 * knowing exact enum type. This class contains int and is convertible to int.
 *
 * Usage example:
 * @CODE
 *
 *      // in some general read function that do not know BCTypeEnum
 *      int bc_type_int = record.val<Enum>("bc_type_selection_key");
 *      ...
 *      // outside of general function
 *      enum { dirichlet, neumann, newton } BCTypeEnum;
 *      BCTypeEnum bc_type = bc_typ_int;
 * @ENDCODE
 *
 */
class Enum {
public:
    Enum() : val_(0) {}
    Enum(int v) :val_(v) {}
    operator int() const {return val_;}
    operator unsigned int() const {return val_;}
private:
    int val_;
};


// Forward declaration
class IteratorBase;
template <class T> class Iterator;

class JSONToStorage;

/**
 * Class for storing and formating input address of an accessor (necessary for input errors detected after readed).
 *
 * To get full path of an accessor we need:
 * - root Input::Type
 * - whole path through the storage
 *
 * TODO:
 * - allow Address with NULL pointers, allow default constructor
 * - How we can get Address with NULL pointer to storage? (currently we need Array::empty_storage_)
 *   - default constructor (should be called only by empty accessors)
 *     see if we can not get empty accessor in json_to_storage
 *     => empty address is error in program
 *   - test NULL pointers in Address::Address(.., ..) constructor
 *   - down ( pokud storage nemuze vracet null , tak zde take nedostaneme null)
 *
 * - find all places where we use Address::storage_head(), check NULL pointer there
 * - where we need StorageArray::new_item() and if we can replace it by add_item()
 */
class Address {
protected:
    struct AddressData {
        /**
         * Pointer to data of parent node in the address tree
         */
        AddressData * parent_;
        /**
         * Index in StorageArray of the parent_ to get actual node.
         */
        unsigned int descendant_order_;
        /**
         * Root of the Input::Type tree.
         */
        const Input::Type::TypeBase *root_type_;
        /**
         * Root of the storage tree.
         */
        const StorageBase *root_storage_;
        /**
         * Actual storage - tip of the storage tree
         */
        const StorageBase * actual_storage_;

        /**
         * Delete whole storage tree when last root input accessor is destroyed.
         */
        ~AddressData() {
        	if (	parent_ == nullptr
        			&& root_storage_ == actual_storage_
        			&& root_type_ ) {
        		delete root_storage_;
        	}
        }
    };

public:
    /**
     * Empty constructor.
     *
     * Constructor should be called only by empty accessors.
     */
    Address();

    /**
     * Basic constructor. We forbids default one since we always need the root input type.
     */
    Address(const StorageBase * storage_root, const Type::TypeBase *type_root);

    /**
     * Copy constructor.
     * TODO: For optimization we can
     * use one vector of storage pointers shared (using shared_ptr) by all accessors along the path.
     */
    Address(const Address &other);

    /**
     * Dive deeper in the storage tree following index @p idx. Assumes that actual node
     * is an StorageArray, has to be asserted.
     */
    const Address * down(unsigned int idx) const;

    /**
     * Getter. Returns actual storage node.
     */
    inline const StorageBase * storage_head() const {
    	ASSERT(data_->actual_storage_, "NULL pointer to storage in address object!!! \n");

    	return data_->actual_storage_;
    }

    /**
     * Produce a full address, i.e. sequence of keys and indices separated by '/',
     * that leads from the root storage and root Input::Type::TypeBase to the actual node in the storage
     * that is nodes_[actual_node_].
     */
    std::string make_full_address() const;

protected:

    /**
     * Shared part of address.
     */
    boost::shared_ptr<AddressData> data_;


};



/**
 * Address output operator.
 */
inline std::ostream& operator<<(std::ostream& stream, const Address & address) {
	return stream << address.make_full_address();
}


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
	typedef ::Input::Type::Record InputType;
    /**
     * Default constructor.
     *
     * Constructor uses empty Address which causes error in program, Address has to be filled.
     */
    Record();

    /**
     * Copy constructor.
     */
    Record(const Record &rec);

    /**
     * Constructs the accessor providing pointer \p store to storage node with list of data of the record and
     * type specification of the record given by parameter \p type.
     */
    Record(const Address &address, const Type::Record type);

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
     * This has similar function as the previous method, but simpler usage in some cases. You has to provide reference to the variable @p value
     * where the value of an optional @p key should be placed. If the key in not present in the input the value of @p value is not changed
     * and the method returns false. If the key has a value the method returns true. Typical usage:
     * @code
     * double param;
     * string other_param;
     * if (rec.opt_val("optional_param", param) ) {
     *      use_param(param);
     * } else if (rec.opt_val("other_param", other_param) ) {
     *      use_other_param(other_param);
     * } else {
     *      ... error, no value for param
     * }
     * @endcode
     */
    template <class Ret>
    inline bool opt_val(const string &key, Ret &value) const;

    /**
     * Returns true if the accessor is empty (after default constructor).
     * TODO: have something similar for other accessors.
     */
    inline bool is_empty() const
    { return (address_.storage_head() == Address().storage_head()); }

    /**
     * Returns address error info.
     */
    EI_Address ei_address() const ;

    /**
     * Get address as string.
     */
    string address_string() const;



protected:
    /**
     * Set address (currently necessary for creating root accessor)
     */
    void set_address(const Address &address);
    friend class JSONToStorage;

    /// Corresponding Type::Record object.
    Input::Type::Record record_type_ ;

    /// Contains address and relationships with record ancestor
    Address address_;
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
	typedef ::Input::Type::AbstractRecord InputType;

    /**
     * Default constructor creates an empty accessor.
     *
     * Constructor uses empty Address which causes error in program, Address has to be filled.
     */
    AbstractRecord();

    /**
     * Copy constructor.
     */
    AbstractRecord(const AbstractRecord &rec);

    /**
     * Constructs the accessor providing pointer \p store to storage node with list of data of the record and
     * type specification of the record given by parameter \p type.
     */
    AbstractRecord(const Address &address, const Type::AbstractRecord type);

    /**
     * Implicit conversion to the \p Input::Record accessor. You can use \p Input::AbstractRecord in the same
     * way as the \p Input::Record.
     */
    operator Record() const;

    /**
     * Returns particular type selected from input. You can use it to construct particular type.
     *
     * @code
     * class MyClass {
     *      MyClass( Input::Record );
     * }
     *
     * if (abstract_record.type() == MyClass.get_input_type())
     *      my_class = new MyClass(abstract_record);        // here the implicit conversion to Input::Record is used
     * @endcode
     */
    Input::Type::Record type() const;

    /**
     * Returns address error info.
     */
    EI_Address ei_address() const ;

    /**
     * Get address as string.
     */
    string address_string() const;


    /**
     * Construct classes given by TYPE key of AbstractRecord.
     *
     * Method uses Input::Factory class. All constructed classes (representing by descendants
     * of AbstractRecord) must be registered to factory (see Input::Factory class) and must have
     * constructors with same parameters (given by Arguments).
     */
    template<class Type, class... Arguments>
    const std::shared_ptr<Type> factory(Arguments... arguments) const;


private:
    /// Corresponding Type::AbstractRecord object.
    Input::Type::AbstractRecord record_type_ ;

    /// Contains address and relationships with abstract record ancestor
    Address address_;
};



/**
 * @brief Accessor to input data conforming to declared Array.
 *
 * There are two possible ways how to retrieve data from Array accessor. First, you can use generic
 * @p copy_to function to copy the data into a given container. Second, you can get an Iterator<Type>
 * and iterate through the Array. Unfortunately, you have to provide Type to the begin() method so this
 * implementation is not fully compliant with standard library. The reason is that in order to speed up compilation of many
 * classes using input accessors we wouldn't have Input::Array a class template that it can be compiled only once.
 * By this reason one can not use BOOST_FOREACH to iterate over Input::Array.
 * TODO: Make Input::Array<Type> wrapper which is compliant with standard library.
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

	typedef ::Input::Type::Array InputType;

    /**
     * Default constructor, empty accessor.
     *
     * Constructor uses empty Address which causes error in program, Address has to be filled.
     */
    Array();

    /**
     * Copy constructor.
     */
    Array(const Array &ar);

    /**
     * Constructs the accessor providing pointer \p store to storage node with list of data of the record and
     * type specification of the record given by parameter \p type.
     */
    Array(const Address &address, const Type::Array type);

   /**
    * Returns iterator to the first element of input array. The template parameter is C++ type you want to
    * read from the array. Only types supported by Input::Interface::Iterator can be used.
    */
   template <class ValueType>
   inline Iterator<ValueType> begin() const;

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

   /**
    * Returns true if the accessor is empty (after default constructor).
    * TODO: have something similar for other accessors.
    */
   inline bool is_empty() const
   { return (address_.storage_head() == Address().storage_head()); }

   /**
    * Returns address error info.
    */
   EI_Address ei_address() const ;

   /**
    * Get address as string.
    */
   string address_string() const;

   /// Need persisting empty instance of StorageArray that can be used to create an empty Address.
   static StorageArray empty_storage_;

private:
    /// Corresponding Type::Array.
    Input::Type::Array array_type_ ;

    /// Contains address and relationships with array ancestor
    Address address_;


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
     * Constructor. Creates iterator effectively pointing to data address_->get_storage()->get_item(index),
     * that is parameter @p address points to StorageArray and parameter @p index gives index into this array.
     *
     */
    IteratorBase(const Address &address, const unsigned int index)
    : address_(address), index_(index)
    {}

    /// Comparison of two Iterators. Do no compare types only position in the storage
    inline bool operator == (const IteratorBase &that) const;

    inline bool operator != (const IteratorBase &that) const;

    /**
     * Implicit conversion to bool. Returns true if iterator points to non-null storage.
     */
    inline operator bool() const;

    /**
     * Return index in an array or record.
     */
    inline unsigned int idx() const;

    /**
     * Returns address
     */
    const Address &get_address() const
    { return address_; }


protected:
    Address address_;
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
     * Iterator is not default constructible.
     *
     * Constructor uses empty Address which causes error in program, Address has to be filled.
     */
    Iterator() : IteratorBase( Address(), 0) {}

    /**
     * Constructor. Creates iterator effectively pointing to data address_->get_storage()->get_item(index),
     * that is parameter @p address points to StorageArray and parameter @p index gives index into this array.
     * Parameter @p type is Input::Type of object the iterator points to.
     *
     *
     */
    Iterator(const Input::Type::TypeBase &type, const Address &address, const unsigned int index)
    : IteratorBase(address, index), type_( type_check_and_convert(type))
    {}

    /// Prefix. Advance operator.
    inline Iterator<T> &operator ++ ();

    /// Prefix. Back operator.
    inline Iterator<T> &operator -- ();

    /**
     *  Dereference operator * ; Shouldn't we return type T, i.e. try to cast from OutputType to T ??
     */
    inline OutputType operator *() const;

    /**
     *  Dereference operator can be used only for iterators to accessors Record, AbstractRecord, and Array.
     */
    inline OutputType *operator ->() const;


private:

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
    static inline ReadType value(const Address &a, const InputType&) { return ReadType( a.storage_head()->get_int() ); }
};


template<>
struct TypeDispatch<int> {
    typedef Input::Type::Integer InputType;
    typedef const int ReadType;
    typedef int TmpType;
    static inline ReadType value(const Address &a, const InputType&) { return a.storage_head()->get_int(); }
};

template<>
struct TypeDispatch<bool> {
    typedef Input::Type::Bool InputType;
    typedef const bool ReadType;
    typedef int TmpType;
    static inline ReadType value(const Address &a, const InputType&) { return a.storage_head()->get_bool(); }
};

template<>
struct TypeDispatch<double> {
    typedef Input::Type::Double InputType;
    typedef const double ReadType;
    typedef int TmpType;
    static inline ReadType value(const Address &a, const InputType&) { return a.storage_head()->get_double(); }
};


template<>
struct TypeDispatch<string> {
    typedef Input::Type::String InputType;
    typedef const string ReadType;
    typedef int TmpType;
    static inline ReadType value(const Address &a, const InputType&) { return a.storage_head()->get_string(); }
};


template<>
struct TypeDispatch<AbstractRecord> {
    typedef Input::Type::AbstractRecord InputType;
    typedef AbstractRecord ReadType;
    typedef AbstractRecord TmpType;
    static inline ReadType value(const Address &a, const InputType& t) { return AbstractRecord(a, t); }
};


template<>
struct TypeDispatch<Record> {
    typedef Input::Type::Record InputType;
    typedef Record ReadType;
    typedef Record TmpType;
    static inline ReadType value(const Address &a, const InputType& t) { return Record(a,t); }
};


template<>
struct TypeDispatch<Array> {
    typedef Input::Type::Array InputType;
    typedef Array ReadType;
    typedef Array TmpType;
    static inline ReadType value(const Address &a, const InputType& t) { return Array(a,t); }

};

template<>
struct TypeDispatch<FilePath> {
    typedef Input::Type::FileName InputType;
    typedef FilePath ReadType;
    typedef int TmpType;
    static inline ReadType value(const Address &a, const InputType& t) { return FilePath(a.storage_head()->get_string(), t.get_file_type() ); }

};




} // closing namespace internal





} // closing namespace Input


// include implementation of templates and inline methods
#include "accessors_impl.hh"



#endif
