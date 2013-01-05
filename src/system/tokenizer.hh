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
 * reading of the text files. Actual tokenizer use backslash '\\' as the escape character, double quotas '"'as quotation
 * character, and space ' ' or tabelator '\\t' as the separator of tokens.
 *
 * !! Used token separator @p escaped_list_separator do not provide possibility to merge several consecutive
 * separator characters into one separator. Consequently, there appears empty tokens when there more spaces
 * then one separating tokens. To overcome this, we drop every empty token.
 *	
 * Provides:
 * - method to read @p next_line, automatically skipping empty lines
 * - iterating over tokens on current line
 * - number of lines that the tokenizer has read -- method line_num
 *
 * Example of usage:
 * @code
 * Tokenizer(in_stream);
 * @endcode
 *
 * TODO:
 * - method to reopen - skip to the beginning
 *
 *
 */
class Tokenizer {
public:
    /**
     * Shortcut for boost tokenizer.
     */
    typedef boost::escaped_list_separator<char> Separator;
    //typedef boost::tokenizer<boost::char_separator<char> > BT;
    typedef boost::tokenizer<Separator> BT;

    /**
     * Opens a file given by file path @p fp. And construct the tokenizer over the
     * input stream for this file.
     * The stream is read from its actual position. The separator of the tokens is
     * either tabelator '\\t' or space ' '.
     *
     */
    Tokenizer(const  FilePath &fp);
    /**
     * Construct the tokenizer over given input stream @p in.
     * The stream is read from its actual position. The separator of the tokens is
     * either tabelator '\\t' or space ' '.
     */
    Tokenizer( std::istream &in);
    /**
     * Skip forward to the line that match given string. 
     * The tokenizer is set to the begining of that line.
     * Returns true if the @p pattern has been found before the end of file.
     * Optionally, if the parameter @p end_search_pattern is provided, the search is
     * stopped after search for @p pattern on the line where @p end_search_pattern was detected.
     * Next line is not read.
     * 
     * TODO: similar method that use regular expressions (e.g. from boost)
     * TODO: add option to find the pattern in the whole file (wrap over the end of file)
     */
    bool skip_to(const std::string &pattern, const std::string &end_search_pattern="");
    
    /**
     * Drops remaining tokens on the current line and reads the new one.
     * A warning is reported in the case of unprocessed tokens.
     * The lines without any tokens are skipped, but counted into
     * number reported by @p line_num. Retuns false if we reach the end of file
     * otherwise returns true.
     *
     * Optional parameter @p assert_for_remaining_tokens can be set false if you
     * want to ignore remaining tokens on current line. Otherwise an warning for the user is
     * produced since possibly there is error in the data format.
     */
    bool next_line(bool assert_for_remaining_tokens=true);
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
      // skip empty tokens (consecutive separators)
      while (! eol() && (*tok_).size()==0 ) {position_++; ++tok_;}
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
     * Returns position on line.
     */
    inline unsigned int pos() const
        { return position_;}

    /**
     * Returns number of lines read by the tokenizer.
     * After first call of @p next_line this returns '1'.
     */
    inline unsigned int line_num() const
        {return line_counter_;}

    /**
     * Returns file name.
     */
    inline const std::string &f_name() const
        {return f_name_;}

    /**
     * Returns full position description.
     */
    std::string position_msg() const;

    /**
     * Destructor close the file if it was opened by tokenizer itself.
     */
    ~Tokenizer();

private:
    // reset tokenizer for actual line
    void set_tokenizer();
    
    /// File name (for better error messages)
    std::string f_name_;
    /// Pointer to internal stream , if tokenizer is constructed form FilePath object.
    std::ifstream *own_stream_;
    /// Input stream.
    std::istream *in_;
    /// Current line
    std::string line_;
    /// Number of liner read by the tokenizer.
    unsigned int line_counter_;
    unsigned int position_;

    /// Line token iterator
    BT::iterator tok_;
    /// Separator function used by the tokenizer
    Separator separator_;
    /// Line tokenizer (container like object).
    BT line_tokenizer_;
};




#endif /* TOKENIZER_HH_ */
