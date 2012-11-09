/*
 * tokenizer.hh
 *
 *  Created on: Nov 9, 2012
 *      Author: jb
 */

#ifndef TOKENIZER_HH_
#define TOKENIZER_HH_

#include <istream>
#include <boost/tokenizer.hpp>

using namespace std;


class Tokenizer {
public:
    typedef boost::tokenizer<boost::char_separator<char> > BT;

    Tokenizer( istream &in);
    void next_line();
    inline const std::string & operator *() const
        { return *tok_; }
    inline BT::iterator & operator ++()
        {position++; return ++tok_;}
    inline unsigned int line_num() const
        {return line_counter_;}
private:
    istream &in_;
    string line_;
    unsigned int line_counter_;

    unsigned int position;
    BT::iterator tok_;
    BT line_tokenizer_;
};




#endif /* TOKENIZER_HH_ */
