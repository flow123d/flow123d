/*
 * file_name.cc
 *
 *  Created on: May 23, 2012
 *      Author: jb
 */

#include "file_path.hh"
#include "system/system.hh"
#include "system/exceptions.hh"


// static data members
map<string,string> FilePath::placeholder;
string FilePath::output_dir="";
string FilePath::root_dir="";

FilePath::FilePath(const string file_path, const  FileType ft) {
    if (output_dir == "") {
        xprintf(Warn, "Creating FileName object before set_io_dirs is called. No file path resolving.\n");
        abs_file_path = file_path;
        return;
    }

    if (ft == input_file) {
        abs_file_path = root_dir + DIR_DELIMITER + file_path;
        substitute_value();
    } else if (ft == output_file) {
        if (file_path[0] == DIR_DELIMITER) {
            THROW( ExcAbsOutputPath() << EI_Path( file_path ) );
        }
        abs_file_path = output_dir + DIR_DELIMITER + file_path;
        substitute_value();
    }
}



void FilePath::set_io_dirs(const string working_dir, const string root_input_dir,const string input,const string output) {
    // root directory
    root_dir = root_input_dir;

    // relative output dir is relative to working directory
    // this is possibly independent of position of the main input file
    if (output[0] == DIR_DELIMITER) output_dir = output;
    else output_dir = working_dir + DIR_DELIMITER + output;

    // the relative input is relative to the directory of the main input file
    add_placeholder("${INPUT}", input);
}



void FilePath::add_placeholder(string key,string value) {
    placeholder[key] = value;
}


void FilePath::substitute_value() {
    for (std::map<std::string,std::string>::const_iterator it = this->placeholder.begin(); it != this->placeholder.end(); ++it) {
        size_t i = abs_file_path.find(it->first,0);
        if (i != std::string::npos) {
            if (it->second == "" ) xprintf(Warn, "Substituting placeholder %s with empty value.\n", it->first.c_str());
            abs_file_path.replace(i, it->first.size(), it->second);
        }
    }
}
