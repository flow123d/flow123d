/**
 * @defgroup input_types  Input Types 
 *
 * The purpose of the classes in the \p Input::Type namespace is to declare structure of the 
 * input data. There are few basic data types (scalar types), further there is Array type
 * to declare sequences of data of common type. Similarly to C++ struct type, we provide 
 * type Record which is set of data of possibly different types indexed by string keys.
 * Finally, we introduce kind of polymorphism through \p AbstractRecord type that mimics
 * abstract classes. 
 * 
 * Types that are simple to construct (\p Array and scalar types) are 
 * cloned when the copy constructor is called (e.g. in Record::declare_key) On the other hand,
 * some other types have nontrivial initialization (Record, AbstractRecord, Selection).
 * The latter group is not cloned on copy operations, since these are mere
 * \p boost::shared_ptr to structure with actual type. The copy constructor makes only copy of the 
 * shared pointer. These types are also named and their description can be provided.
 *
 * \section scalar_types Basic scalar types
 *
 * The basic scalar types are \p String, \p Bool, \p Integer, \p Double, and \p FileName. First four
 * types directly corresponds to C++ types. Individual instances of \p String and \p Bool 
 * classes are identical.
 * On the other hand, instances of Integer and Double classes
 * can differ since you can specify the interval of the valid values. 
 * Finally, the type \p FileName is used for initialization of \p FilePath objects and
 * thus can be either for input or for output file.
 *
 * Nontrivial scalar type is \p Selection which mimics and should be used to initialize 
 * C++ enum variables. The Selection object identifies possible integer values (that should 
 * corrspond to enum values) strings that should correspond to
 * names of the values (should be same as names in enum). 
 * 
 * Unfortunately there is no way to get names of an enum type as strings so this information
 * has to be provided through method @p add_value.  
 * Every @p Selection object has particular name and description.
 * Description of individual values is optional.
 *
 * \section array_type Input::Type::Array 
 * The particular \p Array type is given by the type of its elementsand possibly 
 * by the limits for the array size.  All the elements of the array must have same type, but
 * you can use polymorphism provided by AbstractRecord.
 * 
 * 
 * \section record_type Input::Type::Record 
 *
 * One instance of @p Record is like one class definition in C++. 
 * You specify its members calling the method @p declare_key.
 * Key names should be valid C++ keywords. Every key represents a value of arbitrary 
 * but given type. For every key, you can provide a default value and key description.
 *
 *
 * \section abstract_record \ref Input::Type::AbstractRecord
 *
 * Mimics polymorphism of C++ classes. Any Record can be derived (using \p Recrod::derive method)
 * from an AbstractRecord. This has two effects. First, all keys of the AbstractRecord 
 * are also members of the inheriting Record. Second, the actual input data for an AbstractRecord
 * has to specify value of the special key 'TYPE' that provides name of one of descendants of 
 * the AbstractRecord. This mechanism allows to influence type of the input data by the input itself,
 * but only from the declared set of descendants of the AbstractRecord.
 * 
 * 
 * 
 * @ingroup input
 *
 */
