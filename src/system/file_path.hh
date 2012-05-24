/*
 * file_name.hh
 *
 *  Created on: May 23, 2012
 *      Author: jb
 */

#ifndef FILE_NAME_HH_
#define FILE_NAME_HH_

#include <string>

#include "exceptions.hh"

//#include "input/json_to_storage.hh"

using namespace std;




/**
 * Dedicated class for storing path to files. It performs some modification of path to input and output files.
 *
 *
 * IONameHandler stores items formed by the combination of a key value and a mapped value.
 * The pairs are used to substitution in paths of the input and output files.
 *
 * If program is running with "-i" switch the "${INPUT}" variable in the given filename is substituted by the value immediately following this switch.
 * If program is running with "-o" switch all output file paths are prefixed by the value immediately following this switch (instead of root directory flow.ini file)
 *
 *
 * @par Another example usage:
 * @code
 * IONameHandler &io_name_handler = *(IONameHandler::get_instance());
 * io_name_handler.get_input_file_name("relative/path/with/${INPUT}/var/mesh.msh");
 * @endcode
 */

class FilePath {
public:

    //TYPEDEF_ERR_INFO( EI_ErrorAddress, Input::JSONPath);
    TYPEDEF_ERR_INFO( EI_Path, string);
    DECLARE_EXCEPTION( ExcAbsOutputPath, << "Can not set absolute path " << EI_Path::qval << "at address: ??"  );

    enum FileType {
        input_file,
        output_file
    };



    /**
     * Translates the given absolute or relative path to a file @p file_path.
     * - For relative input path prepend directory of main input file.
     * - For every input file substitutes ${INPUT} place holder by the input directory given at command line.
     * - Forbids absolute output paths.
     * - For relative output path prepends it by output direcotry given at the command line.
     */
    FilePath(const string file_path, const  FileType ft);

    /**
     * Set current, input and output dir.
     */
    static void set_io_dirs(const string current,const string input,const string output);

    /**
     * This class is implicitly convertible to string.
     */
    inline operator string() const
        {return abs_file_path;}

    /*!
     * @brief Add new item to place holder.
     *
     * Placeholder is extended by adding a single new item. The item can be used in the name of the input or output file name.
     *
     * @par Example usage:
     * @code
     *      FilePath::add_placeholder_item("${SUBST_VAL}", "path/value");
     * @endcode
     *
     * @param[in]   key Key of new item.
     * @param[in]   val Value of new item.
     */
    static void add_placeholder(string key,string value);


private:
    void substitute_value();


    string abs_file_path;

    /// dictionary of placeholders
    static std::map<string,string> placeholder;

    /// Prefix path for output files.
    static string output_dir;
    /// Prefix path for input files (directory of the main input file).
    static string root_dir;
};

#endif /* FILE_NAME_HH_ */
