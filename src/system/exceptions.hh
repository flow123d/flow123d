/*
 * exceptions.hh
 *
 *  Created on: Apr 10, 2012
 *      Author: jb
 */


#ifndef EXCEPTIONS_HH_
#define EXCEPTIONS_HH_

/**
 * @file Basic exceptions used in Flow123d.
 *
 * We are using boost::exceptions .
 */

#include <boost/exception/all.hpp>

/**
 * Basic exception for Flow123d's exceptions. When deriving particular exceptions always use virtual inheritance:
 * @code
 *      struct my_exception : virtual flow_excepiton {};
 * @endcode
 *
 * TODO: Implement rich error info.
 */
struct ExceptionBase : virtual std::exception, virtual boost::exception
{
    void print_exc_data(std::ostream &out) const;
    virtual void print_info(std::ostream &out) const;
    virtual const char * what () const throw ();
};

/**
 * This should be used as a base class for all exceptions that are due to incorrect input form the program user.
 *
 * TODO: Implement kind error messages for this case.
 */
struct InputException : virtual ExceptionBase {};


#define THROW(whole_exception_expr) \
    BOOST_THROW_EXCEPTION( whole_exception_expr)


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
 *
 * DeclException( SomeFlowException , << "You can provide " << EI_VAL(ErrorCode) << " here.");
 */
#define TYPEDEF_ERR_INFO(EI_Type, Type)                                                 \
namespace internal {struct EI_Type##_TAG;}                                              \
internal::ExcStream & EI_Type##_val(internal::ExcStream & es) {                         \
    es.stream_ << internal::NullOutputEnvelope<Type>                                    \
        ( boost::get_error_info< boost::error_info<internal::EI_Type##_TAG, Type > >(es.exc_), false );     \
    return es;                                                                          \
}                                                                                       \
internal::ExcStream & EI_Type##_qval(internal::ExcStream & es) {                         \
    es.stream_  << internal::NullOutputEnvelope<Type>                                    \
        ( boost::get_error_info< boost::error_info<internal::EI_Type##_TAG, Type > >(es.exc_), true );     \
    return es;                                                                          \
}                                                                                       \
\
typedef boost::error_info< internal::EI_Type##_TAG, Type > EI_Type




/**
 * Macro for declaration of exceptions.
 * Example:
 * @code
 *      TYPEDEF_ERR_INFO( EI_Dim1Mismatch, int)
 *      TYPEDEF_ERR_INFO( EI_Dim2Mismatch, int)
 *      DeclareException( ExcDimensionMismatch, << "Dimensions dim1=" << EI_Dim1Missmatch_val << " and dim2=" << EI_Dim2Mismatch_val << " should be same.");
 * @endcode
 */
#define DECLARE_EXCEPTION( ExcName, Format)                                  \
struct ExcName : virtual ExceptionBase {                                    \
     virtual void print_info(std::ostream &out) const {                     \
         using namespace internal;                                          \
         ExcStream estream(out, *this);                                     \
         estream Format ;                                                   \
         out << std::endl;                                                  \
     }                                                                      \
     virtual ~ExcName() throw () {}                                         \
}

/**
 * Same as previous but the exception is derived from InputException.
 * This should be used for all exceptions due to wrong input from user.
 */
#define DECLARE_INPUT_EXCEPTION( ExcName, Format)                             \
struct ExcName : virtual InputException {                                   \
     virtual void print_info(std::ostream &out) const {                     \
         using namespace internal;                                          \
         ExcStream estream(out, *this);                                     \
         estream Format ;                                                   \
         out << std::endl;                                                  \
     }                                                                      \
     virtual ~ExcName() throw () {}                                         \
}




namespace internal {


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




/*
class ExceptionValue {
public:
    ExceptionValue(const FlowException &exc, bool quoted = false)
    : exc_(exc), quoted_(quoted) {}

    template <class Tag, class Type>
    const NullOutputEnvelope<Type> & operator() ()
    {
        return NullOutputEnvelope<Type>( boost::get_error_info<Tag>(exc_), quoted_ );
    }

private:
    const FlowException &exc_;
    bool quoted_;
};
*/

} // namespace internal
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
