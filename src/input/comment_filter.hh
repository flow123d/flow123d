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
 */
struct uncommenting_fsm : public io::finite_state_machine<uncommenting_fsm> {
    BOOST_IOSTREAMS_FSM(uncommenting_fsm) // Define skip and push.
    typedef uncommenting_fsm self;

    static const int no_comment         = initial_state;
    static const int no_comment_bsl     = initial_state + 1;
    static const int in_quote           = initial_state + 2;
    static const int in_quote_bsl       = initial_state + 3;
    static const int in_comment         = initial_state + 4;
    static const int in_comment_bsl     = initial_state + 5;
    static const int in_comment_bsl_eol_n     = initial_state + 6;
    static const int in_comment_bsl_eol_r     = initial_state + 7;

    typedef boost::mpl::vector<
                row<no_comment,   is<'\\'>,     no_comment_bsl, &self::push>,   // push backslash and any following character
                row<no_comment,   is<'#'>,      in_comment,     &self::skip>,
                row<no_comment,   is<'"'>,      in_quote,       &self::push>,
                row<no_comment,   is_any,       no_comment,     &self::push>,

                row<no_comment,   is_any,       no_comment,     &self::push>,   // push backslash and nay following character

                row<in_comment,   is<'\n'>,     no_comment,     &self::push>,   // comment up to the end line, but do not remove end line
                row<in_comment,   is<'\r'>,     no_comment,     &self::push>,
                row<in_comment,   is<'\\'>,     in_comment_bsl, &self::skip>,   // skip, but possibly skip also whitespace and newlines
                row<in_comment,   is_any  ,     in_comment,     &self::skip>,   // skip everything else in comment

                row<in_comment_bsl, is<'\n'>,   in_comment_bsl_eol_n, &self::skip>,
                row<in_comment_bsl, is<'\r'>,   in_comment_bsl_eol_r, &self::skip>,
                row<in_comment_bsl, is_space,   in_comment_bsl, &self::skip>,
                row<in_comment_bsl, is_any,     in_comment, &self::skip>,

                //same as in_comment, but \r do not exit comment
                row<in_comment_bsl_eol_n,   is<'\n'>,     no_comment,     &self::push>,
                row<in_comment_bsl_eol_n,   is<'\r'>,     in_comment,     &self::push>,
                row<in_comment_bsl_eol_n,   is<'\\'>,     in_comment_bsl, &self::skip>,
                row<in_comment_bsl_eol_n,   is_any,       in_comment,     &self::skip>,

                //same as in_comment, but \n do not exit comment
                row<in_comment_bsl_eol_r,   is<'\n'>,     in_comment,     &self::push>,
                row<in_comment_bsl_eol_r,   is<'\r'>,     no_comment,     &self::push>,
                row<in_comment_bsl_eol_r,   is<'\\'>,     in_comment_bsl, &self::skip>,
                row<in_comment_bsl_eol_r,   is_any  ,     in_comment,     &self::skip>,

                row<in_quote,     is<'"'>,      no_comment,     &self::push>,
                row<in_quote,     is<'\\'>,     in_quote_bsl,   &self::push>,
                row<in_quote,     is_any,       in_quote,       &self::push>,

                row<in_quote_bsl, is_any,       in_quote,       &self::push>
            > transition_table;
};

typedef io::finite_state_filter<uncommenting_fsm> uncommenting_filter;


} // namespace Input


#endif /* COMMENT_FILTER_HH_ */
