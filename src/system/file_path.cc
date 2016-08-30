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
#include <fstream>

#include "file_path.hh"
#include "system.hh"
#include <type_traits>



/**
 * Helper method, create path string from vector of strings
 */
string create_path_from_vec(vector<string> sub_paths) {
	ASSERT_GT(sub_paths.size(), 0).error();

	if (sub_paths.size() == 1) {
		return sub_paths[0];
	} else {
		boost::filesystem::path path = sub_paths[0];
		for (unsigned int i=1; i<sub_paths.size(); ++i)
			path = path / sub_paths[i];
		return path.string();
	}
}


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

    bool is_abs_path = boost::filesystem::path( convert_for_check_absolute(file_path) ).is_absolute();
    substitute_value(file_path);
	abs_file_path_ = std::make_shared<boost::filesystem::path>(file_path);
    if (ft == input_file) {
    	if ( is_abs_path ) {
    	} else {
            abs_file_path_ = std::make_shared<boost::filesystem::path>(*root_dir / file_path);
    	}
    } else if (ft == output_file) {
        if ( is_abs_path ) {
            if (abs_file_path_->string().substr(0, output_dir->string().size()) == output_dir->string()) {
            	abs_file_path_ = std::make_shared<boost::filesystem::path>( abs_file_path_->string().substr(output_dir->string().size()+1) );
            } else {
                THROW( ExcAbsOutputPath() << EI_Path( abs_file_path_->string() ) );
            }
        }
        abs_file_path_ = std::make_shared<boost::filesystem::path>(*output_dir / *abs_file_path_);
    }

}


FilePath::FilePath(vector<string> sub_paths, const  FileType ft)
: FilePath(create_path_from_vec(sub_paths), ft) {}


FilePath::FilePath(string file_path)
: FilePath(file_path, FilePath::output_file) {}


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
	if (boost::filesystem::path( convert_for_check_absolute(root) ).is_absolute()) {
		root_dir = std::make_shared<boost::filesystem::path>(root);
	} else {
		boost::filesystem::path abs_root_path = boost::filesystem::canonical( boost::filesystem::current_path() / root );
		root_dir = std::make_shared<boost::filesystem::path>(abs_root_path);
	}

	// set output directory
	// if param output is absolute output_dir = output
	// else output_dir = root / output
	// the resulting relative path is always completed to absulute path
	boost::filesystem::path full_output_path;
    if ( !boost::filesystem::path( convert_for_check_absolute(output) ).is_absolute() ) {
    	if (boost::filesystem::path( convert_for_check_absolute(root) ).is_absolute()) {
    		full_output_path = boost::filesystem::path(root) / output;

    		create_dir(full_output_path);
    	} else {
    		boost::filesystem::path output_path = boost::filesystem::path(root) / output;
    		create_dir( output_path );
    		full_output_path = boost::filesystem::canonical( boost::filesystem::current_path() / output_path);
    	}
    } else {
    	full_output_path = boost::filesystem::path(output);
    	create_dir(full_output_path);
    }
	output_dir = std::make_shared<boost::filesystem::path>( full_output_path );

	// the relative input is relative to the directory of the main input file
    add_placeholder("${INPUT}", input);
}


string FilePath::set_dirs_from_input(const string main_yaml, const string input, const string output) {
	boost::filesystem::path input_path;
	if ( !boost::filesystem::path( convert_for_check_absolute(main_yaml) ).is_absolute() ) {
    	input_path = boost::filesystem::path(".") / main_yaml;
    } else {
    	input_path = boost::filesystem::path(main_yaml);
    }
    if (! boost::filesystem::exists(input_path))
        THROW(ExcFileOpen() << EI_Path(input_path.string()));

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
        create_dir( abs_file_path_->parent_path() );
    }
}


string FilePath::parent_path() const {
	return abs_file_path_->parent_path().string();
}


string FilePath::filename() const {
	return abs_file_path_->filename().string();
}


string FilePath::stem() const {
	return abs_file_path_->stem().string();
}


string FilePath::extension() const {
	return abs_file_path_->extension().string();
}


string FilePath::cut_extension() const {
	boost::filesystem::path path = abs_file_path_->parent_path() / abs_file_path_->stem();
	return path.string();
}



template <class Stream>
void FilePath::open_stream(Stream &stream) const
{
    if ( std::is_same<Stream, ifstream>::value ) ASSERT(file_type_ == FileType::input_file);
    if ( std::is_same<Stream, ofstream>::value ) ASSERT(file_type_ == FileType::output_file);

    if (file_type_ == FileType::input_file)
        stream.open(abs_file_path_->string().c_str(), ios_base::in);
    else
        stream.open(abs_file_path_->string().c_str(), ios_base::out);

    if (! stream.is_open())
        THROW(ExcFileOpen() << EI_Path(abs_file_path_->string()));

}



bool FilePath::exists() const
{
    return boost::filesystem::exists( *(this->abs_file_path_) );
}



template void FilePath::open_stream(ifstream &stream) const;
template void FilePath::open_stream(ofstream &stream) const;
template void FilePath::open_stream( fstream &stream) const;

string FilePath::convert_for_check_absolute(const string path) {
	ASSERT(path.length()).error("Empty path.");

	if (path[0] == '/') {
		return "/" + path;
	} else {
		return path;
	}
}



FilePath::operator string() const {
	return abs_file_path_->string();
}


void FilePath::create_dir(const boost::filesystem::path &dir)
{
    try {
        boost::filesystem::create_directories(dir);
    } catch (boost::filesystem::filesystem_error &e) {
        THROW(ExcMkdirFail() << EI_Path( dir.string() ) << make_nested_message(e) );
    }
}


bool FilePath::operator ==(const FilePath &other) const
    {return *abs_file_path_ == boost::filesystem::path( string(other) ); }


std::ostream& operator<<(std::ostream& stream, const FilePath& fp) {
    return ( stream << string(fp) );
}
