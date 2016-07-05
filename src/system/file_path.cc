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

#include "file_path.hh"
#include "system.hh"


// static data members
map<string,string> FilePath::placeholder;
std::shared_ptr<boost::filesystem::path> FilePath::output_dir=std::make_shared<boost::filesystem::path>(".");
std::shared_ptr<boost::filesystem::path> FilePath::root_dir=std::make_shared<boost::filesystem::path>(".");

FilePath::FilePath(string file_path, const  FileType ft)
: file_type_(ft)
{
    if (*output_dir == boost::filesystem::path(".")) {
    	WarningOut() << "Creating FileName object before set_io_dirs is called. No file path resolving." << std::endl;
    	abs_file_path_ = std::make_shared<boost::filesystem::path>(file_path);
        return;
    }

    substitute_value(file_path);
	abs_file_path_ = std::make_shared<boost::filesystem::path>(file_path);
    if (ft == input_file) {
    	if ( abs_file_path_->is_absolute() ) {
    	} else {
            abs_file_path_ = std::make_shared<boost::filesystem::path>(*root_dir / file_path);
    	}
    } else if (ft == output_file) {
        if ( abs_file_path_->is_absolute() ) {
            if (abs_file_path_->string().substr(0, output_dir->string().size()) == output_dir->string()) {
            	abs_file_path_ = std::make_shared<boost::filesystem::path>( abs_file_path_->string().substr(output_dir->string().size()+1) );
            } else {
                THROW( ExcAbsOutputPath() << EI_Path( abs_file_path_->string() ) );
            }
        }
        abs_file_path_ = std::make_shared<boost::filesystem::path>(*output_dir / *abs_file_path_);
    }

}


FilePath::FilePath()
    : abs_file_path_( std::make_shared<boost::filesystem::path>("/__NO_FILE_NAME_GIVEN__") ),
      file_type_(output_file)
{}



void FilePath::set_io_dirs(const string working_dir, const string root, const string input, const string output) {
	ASSERT_EQ(working_dir, ".").error();
    FilePath::set_dirs(root, input, output);
}


void FilePath::set_dirs(const string root, const string input, const string output) {
    // root directory
	root_dir = std::make_shared<boost::filesystem::path>(root);

	// set output directory
	// if param output is absolute output_dir = output
	// else output_dir = root / output
	// the resulting relative path is always completed to absulute path
	boost::filesystem::path output_path = boost::filesystem::path(output);
    if ( !output_path.is_absolute() ) {
    	boost::filesystem::path root_output_path = boost::filesystem::path(root) / output_path;
    	boost::filesystem::create_directories(root_output_path);
    	if ( !root_output_path.is_absolute() ) {
        	boost::filesystem::path abs_full_path = boost::filesystem::canonical( boost::filesystem::current_path() / root_output_path );
        	output_dir = std::make_shared<boost::filesystem::path>(abs_full_path);
    	} else {
    		output_dir = std::make_shared<boost::filesystem::path>(root_output_path);
    	}
    } else {
    	output_dir = std::make_shared<boost::filesystem::path>(output_path);
    }

    // the relative input is relative to the directory of the main input file
    add_placeholder("${INPUT}", input);
}


string FilePath::set_dirs_from_input(const string main_yaml, const string input, const string output) {
	boost::filesystem::path input_path(main_yaml);
	if ( !input_path.is_absolute() ) {
    	input_path = boost::filesystem::path(".") / main_yaml;
    }
	FilePath::set_dirs(input_path.parent_path().string(), input, output);

    return input_path.filename().string();
}



void FilePath::add_placeholder(string key,string value) {
    placeholder[key] = value;
}


void FilePath::substitute_value(string &path) {
    for (std::map<std::string,std::string>::const_iterator it = this->placeholder.begin(); it != this->placeholder.end(); ++it) {
        size_t i = path.find(it->first,0);
        if (i != std::string::npos) {
            ASSERT(it->second != "")(it->first).warning("Substituting placeholder with empty value.");
            path.replace(i, it->first.size(), it->second);
        }
    }
}


const string FilePath::get_absolute_working_dir() {
	return boost::filesystem::current_path().string() + "/";
}



void FilePath::create_output_dir() {
    if (file_type_ == output_file) {
        boost::filesystem::create_directories(
                abs_file_path_->parent_path()
                );
    }
}


FilePath::operator string() const {
	return abs_file_path_->string();
}


bool FilePath::operator ==(const FilePath &other) const
    {return *abs_file_path_ == boost::filesystem::path( string(other) ); }
