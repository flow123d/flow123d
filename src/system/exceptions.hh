/*
 * exceptions.hh
 *
 *  Created on: Apr 10, 2012
 *      Author: jb
 */


#ifndef EXCEPTIONS_HH_
#define EXCEPTIONS_HH_


#include <boost/exception/all.hpp>
#include <iostream>
#include <string>



/**
 * Macro for throwing with saving place of the throw.
 * Usage:
 * @code
 *      THROW( ExcError() << EI_SomeValue(42) );
 * @endcode
 * EI_SomeValue is a error_info object for transfer of values form throw point to catch point. See EI<Tag,Type> class template.
 */
#define THROW(whole_exception_expr) BOOST_THROW_EXCEPTION( whole_exception_expr)


namespace internal {
    class ExcStream;
}

/**
 * Basic exception class. We use boost::exception as parent in order to allow passing
 * some data through the exception object from the throw point to the catch point.
 * See DECLARE_EXCEPTION macro for usage.
 *
 * When deriving particular exceptions always use virtual inheritance:
 * @code
 *      struct my_exception : virtual flow_excepiton {};
 * @endcode
 *
 */
class ExceptionBase : public virtual std::exception, public virtual boost::exception
{
public:
    /// Default constructor, calls fill_stacktrace.
    ExceptionBase();
    /// Copy constructor, performs deep copy of stacktrace.
    ExceptionBase(const ExceptionBase &other);
    /// Call GNU backtrace if available, save call stack information into @p stacktrace member.
    void fill_stacktrace();
    /// Prints formated stacktrace into given stream @p out.
    void print_stacktrace(std::ostream &out) const;
    /**
     * Purely virtual method, that should be implemented by descendants. Prints specific erro message into
     * stream @p out. In particular you can use macros DECLARE_EXCEPTION or INPUT_EXCEPTION for easy decalrations.
     */
    virtual void print_info(std::ostringstream &out) const=0;
    /**
     *  Overloaded method for output the exception message if it is not catched.
     *  Implements composition of complex message including diagnostic informations and stack trace.
     *  Should not be overloded in descendant classes. Use @p print_info instead.
     */
    virtual const char * what () const throw ();
    /// Destructor, possibly free stacktrace.
    virtual ~ExceptionBase() throw ();

private:

    /// Array of backtrace frames returned by glibc backtrace_symbols.
    char ** stacktrace;

    /// Size of stacktrace table - number of frames.
    int n_stacktrace_frames;
};

/**
 * Base class for "input exceptions" that are exceptions caused by incorrect input form the user
 * not by some internal error.
 */
class InputException : public virtual ExceptionBase
{
public:
    virtual const char * what () const throw ();
    virtual ~InputException() throw ();

};



/**
 * Macro for declaration of exceptions. You can use error info classes (EI<Tag,Type> ) and its modifiers to
 * output relevant data.
 * Example:
 * @code
 *      TYPEDEF_ERR_INFO( EI_Dim1Mismatch, int);
 *      TYPEDEF_ERR_INFO( EI_Dim2Mismatch, int);
 *      DECLARE_EXCEPTION( ExcDimensionMismatch, << "Dimensions dim1=" << EI_Dim1Missmatch::val << " and dim2=" << EI_Dim2Mismatch::val << " should be same.");
 * @endcode
 */
#define DECLARE_EXCEPTION( ExcName, Format)                                 \
struct ExcName : public virtual ::ExceptionBase {                           \
     virtual void print_info(std::ostringstream &out) const {               \
         using namespace internal;                                          \
         ::internal::ExcStream estream(out, *this);                         \
         estream Format ;                                                   \
         out << std::endl;                                                  \
     }                                                                      \
     virtual ~ExcName() throw () {}                                         \
}

/* ExcName const &_exc=*this; */


/**
 * Same as previous but the exception is derived from InputException.
 * This should be used for all exceptions due to wrong input from user.
 */
#define DECLARE_INPUT_EXCEPTION( ExcName, Format)                             \
struct ExcName : public virtual ::InputException {                                   \
     virtual void print_info(std::ostringstream &out) const {                     \
         using namespace internal;                                          \
         ::internal::ExcStream estream(out, *this);                                     \
         estream Format ;                                                   \
         out << std::endl;                                                  \
     }                                                                      \
     virtual ~ExcName() throw () {}                                         \
}



/**
 * @brief Macro to simplify declaration of error_info types.
 *
 * These are  used to pass data through boost exceptions from the throw point to the catch point, or possibly collect
 * various data along stack rewinding when an exception is thrown. It decalres type ErrorCode_EI, that can be used to
 * pass data into an exception and onother type ErrorCode, that can be used to refer these date in formats of error messages in declaration of exceptions.
 * example of usage:
 *
 * Notes (updates):
 * Tag is only for uniqueness, the only important is typedef name.
 * This macro also provides manipulators Tag_val and Tag_qval for simplifying formationg of error messages.
 * When declaring an exception inside a class (may be practical if you want to pass the instance of the class through the exception),
 * you have to add the class name to the EI_... names, e.g.
 * @code
 * class Class {
 *      TYPEDEF_ERR_INFO(EI_Class, Class);                      // passing the Class itself
 *      DECLACRE_EXCEPTION(ExcClass, << Class::EI_Class_val);   //
 * }
 * ostream & operator <<(ostream &, const & Class) ...          // of course you have to provide output operator for the Class
 * @endcode
 *
 * @code
 * TYPEDEF_ERR_INFO( ErrorCode, int) // declares type ErrorCode_EI and ErrorCode
 *
 * ...
 *
 * throw SomeException() << ErrorCode_EI(10); // here you pass 'int'
 *
 * ...
 *
 * catch (SomeException & exception) {
 *      int const * = boost::get_error_info<TAG_ErrorCode_EI>(); // here you get pointer to const 'int'
 * }
 *
 * @endcode
 *
 * or you rather declare an exception
 * @code
 * DeclException( SomeFlowException , << "You can provide " << EI_VAL(ErrorCode) << " here.");
 * @endcode
 * 
 * TODO: static modifier is necessary if we define inside a class but this limits declaration of
 * any exception that use the Tag to the same compilation unit !! Better solution ?
 *
 * TODO: have method that accepts an lambda function that can be used to extract particular value from stored object, something like:
 *
 * TYPEDEF_ER_INFO( EI_Value, type_info * )
 * DECLARE_EXCEPTION( Exc, << EI_Value_val(*this, (_1)->name()) );
 *
 * This way we can print "NO_VALUE" if eoor_info is not passed to the exception, but if the value is given
 * we can apply the lambda expression to extract name of an type_info object.
 *
 *
 * I didn't foud a way how simplify output of more complex expressions like
 *
 * DECLARE_EXCEPTION(Exc, << NullOutputEnvelope( EI_Value::ptr(_exc) ? &( (*EI_Value::ptr(_exc))->name() ) : NULL ) );
 *
 * !!! DO NOT WORK SINCE you can not take address (& operator) of a result of some method.
 * We need temporary variable, guessing that this can not be done without helper function, i.e. lambda function.
 *
 * There are several problems:
 * 1) To get the pointer to the value, we need reference to the exception. In the case of EI<>::val manipulators
 *    this information is provided from the ExcStream, but here the pointer appears inside an expression so we have to
 *    provide the exception (_exc is local variable containing reference to current exception )
 *
 * 2) Say, we have a pointer PTR and EXPR(x) is an expression we want to evaluate. We want to
 *    print out "NO_VALUE" if PTR == NULL or print result of EXPR(x) if it is not.
 *    This can only be done by passing a lambda function. The boost::lambda is not sufficient since
 *    it supports only very limited set of overloaded operators, namely there is no way how to get
 *    members of a class. The only way is  C++11 standard which provides true lambda functions, but there
 *    you have to provide type of parameters and thus we need EI_Value twice, like
 *
 *     EI_Value::expr( _exc, []( EI_Value::type &ref){return ref->name()} )
 *
 *    Which is longer previous expression.
 *    I'm not sure if you need not to define also the resulting type, but I've seen some examples without it,
 *    so possibly compiler can deduce it from parameters. The EI_Value::expr should be something like:
 *
 *    template<class Func>
 *    EI_Value::expr( ExceptionBase &e, Func &f) { return NullOutputEnvelope( ptr(e) ? &( f(*ptr(e)) ) : NULL ); }
 *
 *
 *
 *
 */
#define TYPEDEF_ERR_INFO(EI_Type, Type)       typedef EI< struct EI_Type##_TAG, Type > EI_Type

/**
 * This class should not be used directly but through macro TYPEDEF_ERR_INFO.
 * It is derived from boost::error_info<tag, type> and similarly as its parent it
 * is tailored for passing a value of type @p Type from the throw point to
 * the catch point. Compared to boost::error_info it provides manipulators @p val and @p qval
 * which can by used in formating exception message to the ExcStream. The first manipulator evaluates
 * directly to the output of the stored value, while @p qval puts the output into single quotas.
 *
 * The static function @p ref can be used when you want to extract and output some particular information
 * from the passed value. However, if no value is given at throw point this function simply aborts. There is
 * probably no way how to make the check and still keep flexibility in manipulation with the result of @p ref
 * function.
 *
 * For usage see documentation of @p TYPEDEF_ERR_INFO mecro.
 */
template<class Tag, class Type>
class EI : public boost::error_info< Tag, Type > {
public:
    typedef typename boost::error_info< Tag, Type> ErrorInfo;

    /// Construction from given value, that has to bee passed to the catch point.
    EI(Type const & value) : ErrorInfo(value) {}
    /**
     * Stream manipulator used to output the stored value. We have to use special stream ExcStream, that
     * has overloaded << operator in order to support manipulators 'val' and 'qval'.
     */
    static internal::ExcStream & val(internal::ExcStream & es);
    /**
     * Stream manipulator for output of quoted value.
     */
    static internal::ExcStream & qval(internal::ExcStream & es);

    /**
     * Returns reference to stored value in given exception object @p e.
     * Check validity of the value.
     */
    static Type const & ref( ExceptionBase const &e);

    /**
     * Similar to the previous but returns pointer to the stored value and check nothing.
     */
    static Type const * ptr( ExceptionBase const &e);

};








namespace internal {

/*
 * Helper class template. Together with its redirection operator it either outputs value pointed by pointer given to the constructor or
 * , if the pointer is null, outputs string 'NO_VALUE'. If the optional parameter quoted is true, the value is printed in single quotation marks.
 */
template <class Type>
class NullOutputEnvelope {
public:
    NullOutputEnvelope( const Type * x, bool quoted =false)
    : x_(x), quoted_(quoted) {}
    inline bool is_null() const {return x_ == NULL;}
    inline bool is_quoted() const { return quoted_; }
    inline const Type & value() const {return (*x_);}
private:
    const Type * x_;
    bool quoted_;
};

template<class Type>
std::ostream& operator<<
    (std::ostream& stream, const NullOutputEnvelope<Type> & value)
{
   if (value.is_null()) return (stream << "NO_VALUE");
   else if (value.is_quoted()) return (stream << "'" << value.value() << "'");
   else return (stream << value.value() );
}


/**
 * Helper stream that knows stream and exception for which we format the message.
 *
 * TODO: Make std manipulators work!
 */

class ExcStream : public std::ostream {
public:
    ExcStream(std::ostream & stream,const ExceptionBase &exc) : stream_(stream), exc_(exc) {}
    std::ostream & stream_;
    const ExceptionBase & exc_;

    // treat exception manipulators
    ExcStream & operator<<(ExcStream & (*pf) (ExcStream &) )
    {
        return pf(*this);
    }
/*
    // treat other manipulators ( This doesn't work. Why?? )
    template <class T>
    ExcStream & operator<<(T & (*pf) (T &) )
    {
        pf(stream_);
        return (*this);
    }*/
};


template < class T>
ExcStream & operator<<(ExcStream & estream, const T & x)
{
    estream.stream_ << x;
    return estream;
}

template < class Tag, class Type, class Func>
internal::ExcStream & operator<<(internal::ExcStream & estream, typename EI<Tag, Type>::template lambda<Func> const & lambda_func)
{
    if (EI<Tag,Type>::get_value_ptr()) estream.stream_ << "NO_VALUE";
    else estream.stream_ << lambda_func.func_(* EI<Tag,Type>::get_value_ptr());
    return estream;
}


} // namespace internal


/**
 * Assert exception with an string message.
 */
TYPEDEF_ERR_INFO( EI_Message, std::string);
DECLARE_EXCEPTION( ExcAssertMsg, << "Violated Assert! " << EI_Message::val);





/***********************************************************************
 * Implementation of templated method
 */

template <class Tag, class Type>
internal::ExcStream & EI<Tag, Type>::val(internal::ExcStream & es) {
    es.stream_ << internal::NullOutputEnvelope<Type>
        ( ptr(es.exc_), false );
    return es;
}



template <class Tag, class Type>
internal::ExcStream & EI<Tag, Type>::qval(internal::ExcStream & es) {
    es.stream_  << internal::NullOutputEnvelope<Type>
        ( ptr(es.exc_), true );
    return es;
}


template <class Tag, class Type>
Type const & EI<Tag, Type>::ref( ExceptionBase const &e)
    {
        Type const * val_ptr = boost::get_error_info< ErrorInfo > (e);
        if (! val_ptr) {
            // try to printout unfinished ErrStream
            std::cerr << "------------------------------------------------------------------------------\n";
            std::cerr << " Fatal Error - dereferencing null pointer when formating an exception message.\n";
            std::cerr << "------------------------------------------------------------------------------\n";
            std::cerr << "** Diagnosting Informations **\n";
            std::cerr <<  boost::diagnostic_information_what( e );
            abort();
        }
        return *val_ptr;
    }


template <class Tag, class Type>
Type const * EI<Tag, Type>::ptr( ExceptionBase const &e)
   { return boost::get_error_info< ErrorInfo > (e); }


/**
 * Assertion that results in an exception.
 */
/*
#ifdef DEBUG_ASSERTS

#define ASSERT_THROW( condition, exception ) \
    (condition) ? : throw

#else

#define ASSERT_THROW( condition, exception )
#endif
*/




/*******************************************************************

Myslenky ohledne exceptions, chybovych hlasek a logovani (co dnes zajistuje xprintf

1) Nutno oddelit nechybove logovani (vcetne warningu)  a chybove zpravy (ty pokud je nekdo neodchyti tak vedou na terminate() )

2) Logovani  ted neresim.

3) Jakakoliv chyba by mela radeji vyvolavat vyjimku (lze odchytit, lze pridavat inforamace behem stack unfolding. )
   A to vcetne ASSERT.

4) Neprijemne je, ze exception musim definovat jako zvlastni tridy, zejmena pokud chci mit nejake specialni formaty chybovych hlasek.
   ale navrhl jsem makro a navazujici implementacni sablony, ktere by to mely zjednodusit. (inspirovano deal.ii)

5) bylo by dobre mit par preddefinovanych vyjimek a preddefinovanych assertions.




*/

#endif /* EXCEPTIONS_HH_ */
