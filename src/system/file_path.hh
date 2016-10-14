/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    file_path.hh
 * @brief   
 */

#ifndef FILE_NAME_HH_
#define FILE_NAME_HH_

#include <string>

#include "system/exceptions.hh"

using namespace std;


namespace boost {
    namespace filesystem  {
        class path;
}}




/**
 * @brief Dedicated class for storing path to input and output files.
 *
 * FilePath objects are constructed from given absolute or relative path to the file and its type (input or output).
 * Before you create any instance of the class you have to call static method @p set_io_dirs to set:
 * - working directory of the program (when it was started)
 * - root input directory, i.e. directory of the main input file (given by -s parameter)
 * - input directory, used to replace ${INPUT} placeholder (given by -i parameter)
 * - output directory, where all output files should be placed (given by -o parameter)
 *
 */

class FilePath {
public:

    /**
     * Reporting failure when openning a file.
     */
    TYPEDEF_ERR_INFO( EI_Path, string);
    DECLARE_EXCEPTION( ExcFileOpen, << "Can not open file: " << EI_Path::qval );
    DECLARE_EXCEPTION( ExcAbsOutputPath, << "Can not set absolute path " << EI_Path::qval << " for an output file."  );
    DECLARE_EXCEPTION( ExcMkdirFail, << "Can not create directory: " << EI_Path::qval );

    /// Possible types of file.
    enum FileType {
        input_file,
        output_file
    };

    /**
     * Default constructor, necessary when using  Input::Record::opt_val() to initialize a FilePath.
     */
    FilePath();

    /**
     * Translates the given absolute or relative path to a file @p file_path depending on the file type @p ft.
     *
     * For input files:
     * - For relative path prepend absolute path of the directory of the main input file (root directory).
     * - Replace ${INPUT} place holder with the input directory given at command line.
     *
     * For output files:
     * - Forbids absolute output paths.
     * - For relative output path prepends it by the output directory given at the command line.
     */
    FilePath(string file_path, const  FileType ft);

    /// Same as previous, but create path from vector of strings.
    FilePath(vector<string> sub_paths, const  FileType ft);

    /// Same as previous but implicitly use FileType::output_file
    FilePath(string file_path);

    /**
     * @brief Obsolete method for set input and output directories.
     *
     * Ensures consistency of unit tests.
     *
     * Set:
     * - working directory (used only if the output directory is relative)
     * - root directory (of the main input file)
     * - input directory to replace ${INPUT} place holder
     * - output directory used as prefix to the output files (relative output dirs are relative to the working directory)
     */
    static void set_io_dirs(const string working_dir, const string root, const string input, const string output);

    /**
     * @brief Method for set input and output directories.
     *
     * Set:
     * - root directory (of the main input file)
     * - input directory to replace ${INPUT} place holder
     * - output directory used as prefix to the output files (relative output dirs are relative to the working directory)
     */
    static void set_dirs(const string root, const string input, const string output);

    /**
     * @brief Method for set input and output directories.
     *
     * Same as previous, but in first argument accepts full path of yaml file and returns filename of this yaml file.
     *
     * Set:
     * - root directory (of the main yaml input file)
     * - input directory to replace ${INPUT} place holder
     * - output directory used as prefix to the output files (relative output dirs are relative to the working directory)
     */
    static string set_dirs_from_input(const string main_yaml, const string input, const string output);

    /**
     * This class is implicitly convertible to string.
     */
    operator string() const;

    /*!
     * @brief Add new item to place holder.
     *
     * Placeholder is extended by adding a single new item. The item can be used in the name of the input or output file name.
     * Currently, the only supported placeholder is ${INPUT}.
     *
     * @par Example usage:
     * @code
     *      FilePath::add_placeholder_item("${SUBST_VAL}", "path/value");
     * @endcode
     *
     * @param[in]   key Key of new item.
     * @param[in]   value Value of new item.
     */
    static void add_placeholder(string key,string value);

    /**
     * Return absolute path of actual working directory.
     */
    static const string get_absolute_working_dir();

    /// Equality comparison operators for FilePaths.
    bool operator ==(const FilePath &other) const;


    /**
     * For an output filepath, the directory part (up to last separator) is
     * extracted and all subdirectories are created if doesn't exist yet.
     */
    void create_output_dir();

    /**
     * Return path to file.
     */
    string parent_path() const;

    /**
     * Return name of file with extension.
     */
    string filename() const;

    /**
     * Return name of file without extension.
     */
    string stem() const;

    /**
     * Return extension of file.
     */
    string extension() const;

    /**
     * Return path to file with filename without extension.
     */
    string cut_extension() const;


    /**
     * Open stream for this FilePath.
     * Open mode is determined from the FilePath type.
     */
    template <class Stream>
    void open_stream(Stream &stream) const;

    /**
     * Return true if the FilePath is a file.
     */
    bool exists() const;

private:
    /**
     * Create a directory, and check for exceptions.
     */
    static void create_dir(const boost::filesystem::path &dir);

    /**
     * Substitutes placeholders in @p path.
     */
    void substitute_value(string &path);

    /**
     * @brief Prepare path string for check absolute path.
     *
     * Check first char of path string. If it is slash '/', add second slash char. Two slashes
     * at begin is necessary for correct output of boost::filesystem::path.is_absolute() method
     * for detection absolute path in unix format ("/home/x/y/z") under cygwin.
     */
    static string convert_for_check_absolute(string path);


    /// Final absolute path to the file.
    std::shared_ptr<boost::filesystem::path> abs_file_path_;

    /// File type
    FileType file_type_;

    /// dictionary of placeholders
    static std::map<string,string> placeholder;

    /// Prefix path for output files.
    static std::shared_ptr<boost::filesystem::path> output_dir;

    /// Prefix path for input files (directory of the main input file).
    static std::shared_ptr<boost::filesystem::path> root_dir;
};

/**
 * @brief Allow redirect FilePath to stream.
 */
std::ostream& operator<<(std::ostream& stream, const FilePath& fp);


#endif /* FILE_NAME_HH_ */
