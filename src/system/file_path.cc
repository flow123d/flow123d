/*
 * file_name.cc
 *
 *  Created on: May 23, 2012
 *      Author: jb
 */


#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "file_path.hh"
#include "system.hh"


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
    if (output[0] == DIR_DELIMITER) {
    	vector<string> dirs;
    	boost::split(dirs, output,  boost::is_any_of("/"));
    	output_dir = "";
    	for (vector<string>::iterator it = dirs.begin(); it != dirs.end(); ++it) {
    	    if ( !(*it).size() ) continue;
        	output_dir = output_dir + DIR_DELIMITER + *it;
            if (!boost::filesystem::is_directory(output_dir)) {
            	boost::filesystem::create_directory(output_dir);
            }
        }
    } else {
    	vector<string> dirs;
    	string full_output = working_dir + DIR_DELIMITER + output;
    	boost::split(dirs, full_output, boost::is_any_of("/"));
    	output_dir = "";
    	for (vector<string>::iterator it = dirs.begin(); it != dirs.end(); ++it) {
    	    if ( !(*it).size() ) continue;
    	    if ( !output_dir.size() ) output_dir = *it;
    	    else output_dir = output_dir + DIR_DELIMITER + *it;
            if (!boost::filesystem::is_directory(output_dir)) {
            	boost::filesystem::create_directory(output_dir);
            }
        }

        boost::filesystem::path working_path = boost::filesystem::path(working_dir);
        boost::filesystem::path output_path = boost::filesystem::path(output);

        if (working_dir[0] != DIR_DELIMITER)
        {
        	boost::filesystem::path curr = boost::filesystem::current_path();
        	working_path = boost::filesystem::canonical( curr / working_path );
        }

        boost::filesystem::path full_path = boost::filesystem::canonical( working_path / output_path );

        boost::filesystem::path curr = boost::filesystem::current_path();
    	output_dir = full_path.string();
#ifdef BOOST_WINDOWS_API
        boost::replace_all(output_dir, "\\", "/");
#endif
    }

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
