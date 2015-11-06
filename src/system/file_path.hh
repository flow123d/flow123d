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


    TYPEDEF_ERR_INFO( EI_Path, string);
    DECLARE_EXCEPTION( ExcAbsOutputPath, << "Can not set absolute path " << EI_Path::qval << " for an output file."  );

    /// Possible types of file.
    enum FileType {
        input_file,
        output_file
    };

    /**
     * Default constructor, necessary when using  Input::Record::opt_val() to initialize a FilePath.
     */
    FilePath()
        : abs_file_path_("/__NO_FILE_NAME_GIVEN__"),
          file_type_(output_file)
    {}

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

    /**
     * Set:
     * - working directory (used only if the output directory is relative)
     * - root directory (of the main input file)
     * - input directory to replace ${INPUT} place holder
     * - output directory used as prefix to the output files (relative output dirs are relative to the working directory)
     */
    static void set_io_dirs(const string working_dir, const string root_input,const string input,const string output);

    /**
     * This class is implicitly convertible to string.
     */
    inline operator string() const
        {return abs_file_path_;}

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

    /// Equality comparison operators for regions.
    inline bool operator ==(const FilePath &other) const
        {return abs_file_path_ == string(other); }


    /**
     * For an output filepath, the directory part (up to last separator) is
     * extracted and all subdirectories are created if doesn't exist yet.
     */
    void create_output_dir();

private:
    /**
     * Substitutes placeholders in @p abs_file_path_.
     */
    void substitute_value();


    /**
     * Test if get path is absolute for used operating system.
     */
    static bool is_absolute_path(const string path);


    /**
     * Check if directory stored in output_dir doesn't exist and create its
     */
    static void create_dir(string dir);


    /**
     * Create canonical path of output directory given by relative path.
     */
    static void create_canonical_path(const string working_dir, const string output);


    /// Final absolute path to the file.
    string abs_file_path_;

    /// File type
    FileType file_type_;

    /// dictionary of placeholders
    static std::map<string,string> placeholder;

    /// Prefix path for output files.
    static string output_dir;

    /// Prefix path for input files (directory of the main input file).
    static string root_dir;
};

#endif /* FILE_NAME_HH_ */
