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
 * @brief Wrapper for throw. Saves the throwing point.
 *
 * Macro for throwing with saving place of the throw. Just shortcut for BOOST_THROW_EXCEPTION. Creates boost kind of exception, that may
 * accept further information through redirection.
 *
 * Usage:
 * @code
 *      THROW( ExcError() << EI_SomeValue(42) );
 * @endcode
 *
 * EI_SomeValue is an @p error_info object for transfer of values form throw point to catch point. See @p EI<Tag,Type> class template.
 *
 * @ingroup exceptions
 */
#define THROW(whole_exception_expr) BOOST_THROW_EXCEPTION( whole_exception_expr)


namespace internal {
    class ExcStream;
}

/**
 * @brief Base of exceptions used in Flow123d.
 *
 * We use boost::exception as parent in order to allow passing
 * some data through the exception object from the throw point to the catch point.
 * See DECLARE_EXCEPTION macro for usage.
 *
 * When deriving particular exceptions always use virtual inheritance:
 * @code
 *      struct my_exception : virtual ExceptionBase {};
 * @endcode
 *
 * @ingroup exceptions
 */
class ExceptionBase : public virtual std::exception, public virtual boost::exception
{
public:
    /// Default constructor, just calls @p fill_stacktrace().
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
 * @brief Base of exceptions due to user input.
 *
 * Base class for "input exceptions" that are exceptions caused by incorrect input from the user
 * not by an internal error.
 *
 * @ingroup exceptions
 */
class InputException : public virtual ExceptionBase
{
public:
    virtual const char * what () const throw ();
    virtual ~InputException() throw ();

};



/**
 * @brief Macro for simple definition of exceptions.
 *
 * First parameter, @p ExcName, is name of the exception class to define. Macro expands to
 * the class derived from @p ExceptionBase and implements virtual method @p print_info
 * with output given by @p Format, the second parameter of the macro. This method is used
 * by what() method to produce specific part of the error message.
 *
 * You can use error info classes (see @p TYPEDEFERR_INFO) and stream modifiers
 * @p val, @p qval to output relevant data. Modifier @p val outputs plain value,
 * modifier @p qval outputs value in quotas.
 *
 * Example:
 * @code
 * TYPEDEF_ERR_INFO( EI_Dim1Mismatch, int);
 * TYPEDEF_ERR_INFO( EI_Dim2Mismatch, int);
 * DECLARE_EXCEPTION( ExcDimensionMismatch,
 *     << "Dimensions dim1=" << EI_Dim1Missmatch::val
 *     << " and dim2=" << EI_Dim2Mismatch::val << " should be same.");
 * @endcode
 *
 * One can also pass whole classes through the exception as long as they are default constructable and copy constructable.
 *
 * Example:
 * @code
 * class Matrix;
 * TYPEDEF_ERR_INFO( EI_Matrix, Matrix);
 * DECLARE_EXCEPTION( ExcWrongMatrixState,
 *     << "Matrix state: " << EI_Matrix::ptr(*this)->state() << endl
 *     << "Matrix info: " << EI_Matrix::ref(*this).info() );
 * @endcode
 *
 * The example shows two ways how to call methods of the object of class Matrix passed through the exception. One can either get pointer to the object or
 * reference. The  @p ref method checks that the pointer to the object is not NULL (so that the object was actually passed).
 * The @p ptr method do not perform the check. However, when using one of these methods you have to guarantee that the @p error_info object is passed to the exception
 * at every throw point that use that exception. Otherwise you get an error, meaningful in case of the @p ref method, seg. fault for the @p ptr method.
 *
 * @ingroup exceptions
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
 * @brief Macro for simple definition of input exceptions.
 *
 * Works in the same way as @p DECLARE_EXCEPTION, just define class derived from
 * @p InputException. Meant to be used for exceptions due to wrong input from user.
 *
 * @ingroup exceptions
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
 * Is used to declare types of data that can be passed through exceptions from the throw point to the catch point, or possibly collect
 * various data along stack rewinding when an exception is thrown.
 *
 * Typical usage:
 * @code
 *   // declares type EI_ParticularNumber
 *   TYPEDEF_ERR_INFO( EI_ParticularNumber, int);
 *   // declares exception that use it
 *   DECLARE_EXCEPTION(SomeException, << " Particular number: " << EI_ParticularNumber::val );
 *   // ... some code ...
 *
 *   // here you pass particular number important for
 *   // determination of cause of the exception
 *   THROW( SomeException() << EI_ParticularNumber(10) );
 *
 * @endcode
 * 
 * @ingroup exceptions
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


/**
 * Exception thrown in xprintf function.
 */
TYPEDEF_ERR_INFO( EI_XprintfMessage, std::string);
DECLARE_EXCEPTION( ExcXprintfMsg, << EI_XprintfMessage::val << ":\n");





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
