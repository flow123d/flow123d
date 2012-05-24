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
    ASSERT( output_dir != "", "Creating FileName object before set_io_dirs is called.\n");
    if (ft == input_file) {
        abs_file_path = root_dir + "/" + file_path;
        substitute_value();
    } else if (ft == output_file) {
        if (file_path[0] == '/') {
            THROW( ExcAbsOutputPath() << EI_Path( file_path ) );
        }
        abs_file_path = output_dir + "/" + file_path;
        substitute_value();
    }
}



void FilePath::set_io_dirs(const string root,const string input,const string output) {
    root_dir = root;

    if (output[0] == '/') output_dir = output;
    else output_dir = root_dir + '/' + output;

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
