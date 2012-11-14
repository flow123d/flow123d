/*
 * tokenizer.hh
 *
 *  Created on: Nov 9, 2012
 *      Author: jb
 */

#ifndef TOKENIZER_HH_
#define TOKENIZER_HH_

#include <boost/tokenizer.hpp>
#include <istream>


class FilePath;


/**
 * @brief Simple class for parsing text files.
 *
 * Boost library provides nice tokenizer. The string is viewed as a container of tokens and
 * you can iterate over them. This class simplify the usage of the boost's tokenizer and further simplify
 * reading of the text files.
 *	
 * Provides:
 * - method to read @p next_line, automatically skipping empty lines
 * - iterating over tokens on current line
 * - number of lines that the tokenizer has read -- method line_num
 *
 * Example of usage:
 * @CODE
 * Tokenizer(in_stream);
 *
 * @ENDCODE
 *
 * TODO:
 * - incorporate skip_to method
 *
 *
 */
class Tokenizer {
public:
    /**
     * Shortcut for boost tokenizer.
     */
    typedef boost::tokenizer<boost::char_separator<char> > BT;

    /**
     * Opens a file given by file path @p fp. And construct the tokenizer over the
     * input stream for this file.
     * The stream is read from its actual position. The separator of the tokens is
     * either tabelator '\t' or space ' '.
     *
     */
    Tokenizer( FilePath &fp);
    /**
     * Construct the tokenizer over given input stream @p in.
     * The stream is read from its actual position. The separator of the tokens is
     * either tabelator '\t' or space ' '.
     */
    Tokenizer( std::istream &in);
    /**
     * Skip forward to the line that match given string. 
     * The tokenizer is set to the begining of that line.
     * Returns true if the pattern has been found before end of file.
     * 
     * TODO: similar method that use regular expressions (e.g. from boost)
     */
    bool skip_to(const std::string &pattern);
    
    /**
     * Drops remaining tokens on the current line and reads the new one.
     * A warning is reported in the case of unprocessed tokens.
     * The lines without any tokens are skipped, but counted into
     * number reported by @p line_num. Retuns false if we reach the end of file
     * otherwise returns true.
     */
    bool next_line();
    /**
     * Dereference of the tokenizer iterator. Returns reference to the string
     * that contains current token.
     */
    const std::string & operator *() const;

    /**
     * Moves to the next token on the line.
     */
    inline BT::iterator & operator ++() {
      if (! eol()) {position_++; ++tok_;}
      return tok_;
    }

    /**
     * Returns true if the iterator is over the last token on the current line.
     */
    inline bool eol() const
        { return tok_ == line_tokenizer_.end(); }
        
    /**
     *  Returns true if at the end of the input stream.
     */    
    inline bool eof() const
        { return in_->eof(); }
    /**
     * Returns number of lines read by the tokenizer.
     * After first call of @p next_line this returns '1'.
     */
    inline unsigned int line_num() const
        {return line_counter_;}


    ~Tokenizer();

private:
    void set_tokenizer();
    
    /// File name (for better error messages)
    std::string f_name_;
    /// Pointer to internal stream , if tokenizer is constructed form FilePath object.
    std::istream *own_stream_;
    /// Input stream.
    std::istream *in_;
    /// Current line
    std::string line_;
    /// Number of liner read by the tokenizer.
    unsigned int line_counter_;
    unsigned int position_;

    /// Line token iterator
    BT::iterator tok_;
    /// Line tokenizer (container like object).
    BT line_tokenizer_;
};




#endif /* TOKENIZER_HH_ */
