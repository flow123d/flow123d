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
 * @file    file_path.cc
 * @brief   
 */

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "file_path.hh"
#include "system.hh"


// static data members
map<string,string> FilePath::placeholder;
string FilePath::output_dir="";
string FilePath::input_root_dir="";

FilePath::FilePath(string file_path, const  FileType ft)
: file_type_(ft)
{
	abs_file_path_ = file_path;
    if (output_dir == "") {
    	WarningOut() << "Creating FileName object before set_io_dirs is called. No file path resolving." << std::endl;
        return;
    }

    substitute_value();
    if (ft == input_file) {
    	if ( FilePath::is_absolute_path(abs_file_path_) ) {
    	} else {
            abs_file_path_ = input_root_dir + DIR_DELIMITER + abs_file_path_;
    	}
    } else if (ft == output_file) {
        if ( FilePath::is_absolute_path(abs_file_path_) ) {
            if (abs_file_path_.substr(0, output_dir.size()) == output_dir) {
            	abs_file_path_=abs_file_path_.substr(output_dir.size()+1);
            } else {
                THROW( ExcAbsOutputPath() << EI_Path( abs_file_path_ ) );
            }
        }
        abs_file_path_ = output_dir + DIR_DELIMITER + abs_file_path_;
    }

}



string FilePath::set_io_dirs(const string main_input_file, const string input, const string output, const string working_dir) {
    // relative output dir is relative to working directory
    // this is possibly independent of position of the main input file

	// set output directory
	boost::filesystem::path full_output_path;
    if ( FilePath::is_absolute_path(output) ) {
    	full_output_path = boost::filesystem::path(output);
    } else {
    	full_output_path = boost::filesystem::path(working_dir + DIR_DELIMITER + output);
    }
    boost::filesystem::create_directories(full_output_path);
    output_dir = full_output_path.string();

    if ( !FilePath::is_absolute_path(output_dir) ) {
    	boost::filesystem::path output_path = boost::filesystem::path(output_dir);
    	boost::filesystem::path full_path = boost::filesystem::canonical( boost::filesystem::current_path() / output_path );
    	output_dir = full_path.string();
#ifdef FLOW123D_HAVE_CYGWIN
    	boost::replace_all(output_dir, "\\", "/");
#endif // FLOW123D_HAVE_CYGWIN
    }

    // the relative input is relative to the directory of the main input file
    add_placeholder("${INPUT}", input);

    // set input root directory
	boost::filesystem::path input_path(main_input_file);
	input_root_dir = input_path.parent_path().string();
    return input_path.filename().string();
}



void FilePath::add_placeholder(string key,string value) {
    placeholder[key] = value;
}


void FilePath::substitute_value() {
    for (std::map<std::string,std::string>::const_iterator it = this->placeholder.begin(); it != this->placeholder.end(); ++it) {
        size_t i = abs_file_path_.find(it->first,0);
        if (i != std::string::npos) {
            ASSERT(it->second != "")(it->first).warning("Substituting placeholder with empty value.");
            abs_file_path_.replace(i, it->first.size(), it->second);
        }
    }
}


bool FilePath::is_absolute_path(const string path) {
	if (path.size() == 0) xprintf(UsrErr, "Path can't be empty!\n");
#ifdef FLOW123D_HAVE_CYGWIN
	if (path.size() == 1) return false;
	return isalpha(path[0]) && (path[1] == ':');
#else
	return path[0] == DIR_DELIMITER;
#endif // FLOW123D_HAVE_CYGWIN
}


const string FilePath::get_absolute_working_dir() {
    string abs_path = boost::filesystem::current_path().string() + "/";
#ifdef FLOW123D_HAVE_CYGWIN
    boost::replace_all(abs_path, "\\", "/");
#endif // FLOW123D_HAVE_CYGWIN
	return abs_path;
}



void FilePath::create_output_dir() {
    if (file_type_ == output_file) {
        boost::filesystem::create_directories(
                boost::filesystem::path(abs_file_path_).parent_path()
                );
    }
}

