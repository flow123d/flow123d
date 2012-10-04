#ifndef IO_NAMEANDLER_H
#define IO_NAMEANDLER_H

#include <string>
#include <map>

using namespace std;

/**
 * @brief Singleton class used to change given 'variables' in the file path.
 *
 * IONameHandler stores items formed by the combination of a key value and a mapped value.
 * The pairs are used to substitution in paths of the input and output files.
 *
 * If program is running with "-i" switch the "${INPUT}" variable in the given filename is substituted by the value immediately following this switch.
 * If program is running with "-o" switch all output file paths are prefixed by the value immediately following this switch (instead of root directory flow.ini file)
 *
 * @par Example usage:
 * @code
 * IONameHandler::get_instance()->get_input_file_name("relative/path/with/${VAR}/flow.ini");
 * @endcode
 *
 * @par Another example usage:
 * @code
 * IONameHandler &io_name_handler = *(IONameHandler::get_instance());
 * io_name_handler.get_input_file_name("relative/path/with/${INPUT}/var/mesh.msh");
 * @endcode
 */
class IONameHandler {
public:
    /*!
         * @brief Returns instance of IONameHandler.
         *
         * Class IONameHandler is created as a singleton (Singleton design pattern). Static method get_instance() returns pointer to IONameHandler object.
         *
         * @par Example usage:
         * @code
         * IONameHandler& io_name_handler = *(IONameHandler::get_instance());
         * @endcode
         */
        static IONameHandler* get_instance();
        /*!
         * @brief Returns absolute path to given input file.
         *
         * @param[in] file_name Filename relatively to root directory.
         * @return absolute path to input file
         */
        string get_input_file_name(string file_name);
        /*!
         * @brief Returns absolute path to given output file.
         *
         * @param[in]   file_name Filename relatively to output directory.
         * @return              absolute path to output file
         */
        string get_output_file_name(string file_name);
        /*!
        * @brief Returns value of output_dir variable.
        */
        string get_output_dir();
        /*!
         * @brief Add new item to place holder.
         *
         * IONameHandler placeholder is extended by adding a single new item. The item can be used in the name of the input or output file name.
         *
         * @par Example usage:
         * @code
         * IONameHandler::get_instance()->add_placeholder_item("${SUBST_VAL}", "path/value");
         * @endcode
         *
         * @param[in]   key Key of new item.
         * @param[in]   val Value of new item.
         * @return              always true
         */
        bool add_placeholder_item(string key,string value);
        /*
         * @brief Removes the key (and its corresponding value) from place holder.
         *
         * @param[in] key The key that needs to be removed.
         * @return The value to which the key had been mapped in place holder, or empty string if the key did not have a mapping.
         */
//  string remove_placeholder_item(string key);
private:
  static IONameHandler* instance;
  string root_dir;
  string output_dir;
  std::map<string,string> placeholder;

  /*!
   * @brief Private constructor prevents instantiation from other classes
   */
  IONameHandler() {};
  /*!
   * @brief Private copy constructor
   */
  IONameHandler(IONameHandler const&) {};
  /*!
   * @brief Private assignment operator - can never be called
   */
  IONameHandler& operator=(IONameHandler const&) { return (*this);};
  /*!
   * @brief Initialization of root directory - where is located flow.ini file.
   */
  void initialize_root_dir();
  /*!
   * @brief Initialization of output directory - the value of "-o" command line variable if is set, otherwise where is located flow.ini file.
   */
  void initialize_output_dir();
  /*!
   * @brief initialization of placeholder standard pairs - the value of "-i" command line variable for "${INPUT}" if is set.
   */
  void initialize_placeholder();
  /*!
   * @brief Returns value of root_dir variable.
   */
  string get_root_dir();
  /*!
   * @brief Substitute all variable in the given file path.
   */
  string substitute_value(string file);
};

#endif