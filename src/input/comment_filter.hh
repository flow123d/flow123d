/*
 * comment_filter.hh
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

#ifndef COMMENT_FILTER_HH_
#define COMMENT_FILTER_HH_

#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 30
#include <boost/mpl/vector.hpp>
#include "input/finite_state_filter.hpp"

namespace Input {

namespace io = boost::iostreams;

/**
 * This structure implements a finite state machine for filtering comments of the input file.
 * The comment starts with '#' at any place save in quotas and after backslash. The comment ends be fors non-backslashed end of line.
 * Backslashed end of line consists of the backslash, any sequence of whitespaces
 *
 * We assume, that json parser use backslash as an escape character, so we also have to escape any special meaning
 * of the next character, but we keep backslash on the place.
 * 
 * Don't ask me for syntactic bloat. This struct has to define \p transition_table of the FSM. Access to the rules is
 * provided by methods derived from \p io::finite_state_machine. Then the FSM is  used to implement filter. 
 */
struct uncommenting_fsm : public io::finite_state_machine<uncommenting_fsm> {
    BOOST_IOSTREAMS_FSM(uncommenting_fsm) // Define skip and push event handlers
    inline void slash_and_push(char ch)   // Event handler for a slash that do not open an comment.
        { push('/'); push(ch); }
    typedef uncommenting_fsm self;

    /**
     * Declaration of all states of FSM.
     */
    static const int no_comment         = initial_state;
    static const int one_slash          = initial_state + 1;
    static const int one_line_comment   = initial_state + 2;
    static const int multi_line_comment = initial_state + 3;
    static const int star_in_comment    = initial_state + 4;
    static const int no_comment_bsl     = initial_state + 5;
    static const int in_quote           = initial_state + 6;
    static const int in_quote_bsl       = initial_state + 7;

    /**
     * Declaration of rules of the FSM.
     */ 
    typedef boost::mpl::vector<
    // format of the table:
    //            actual state, input character set, next state, output i.e. pass or skip input character
    
                row<no_comment,   is<'/'>,      one_slash,     &self::skip>,
                row<no_comment,   is<'"'>,      in_quote,       &self::push>,
                row<no_comment,   is<'\\'>,     no_comment_bsl, &self::push>,   // push backslash and any following character
                row<no_comment,   is_any,       no_comment,     &self::push>,

                row<one_slash,    is<'/'>,      one_line_comment, &self::skip>,
                row<one_slash,    is<'*'>,      multi_line_comment, &self::skip>,
                row<one_slash,    is_any,       no_comment,     &self::slash_and_push>,

                row<one_line_comment,   is<'\n'>,     no_comment,     &self::push>,   // comment up to the end line, but do not remove EOL
                row<one_line_comment,   is<'\r'>,     no_comment,     &self::push>,
                row<one_line_comment,   is_any  ,     one_line_comment,     &self::skip>,   // skip everything else in comment

                row<multi_line_comment, is<'*'>,      star_in_comment,  &self::skip>,
                row<multi_line_comment, is<'\n'>,     multi_line_comment,  &self::push>,    // preserve line numbers
                row<multi_line_comment, is<'\r'>,     multi_line_comment,  &self::push>,
                row<multi_line_comment, is_any,       multi_line_comment,  &self::skip>,

                row<star_in_comment, is<'/'>,         no_comment, &self::skip>,
	              row<star_in_comment, is<'\n'>,        multi_line_comment, &self::push>,    // preserve line numbers
	              row<star_in_comment, is<'\r'>,        multi_line_comment, &self::push>,
                row<star_in_comment, is_any,          multi_line_comment, &self::skip>,

                row<in_quote,     is<'"'>,      no_comment,     &self::push>,
                row<in_quote,     is<'\\'>,     in_quote_bsl,   &self::push>,
                row<in_quote,     is_any,       in_quote,       &self::push>,

                row<in_quote_bsl, is_any,       in_quote,       &self::push>,

                row<no_comment_bsl, is_any,       no_comment,       &self::push>
            > transition_table;
};

/**
 * Declare an io filter based on FSM for filtering comments.
 */
typedef io::finite_state_filter<uncommenting_fsm> uncommenting_filter;


} // namespace Input


#endif /* COMMENT_FILTER_HH_ */
