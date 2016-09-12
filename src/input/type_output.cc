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
 * @file    type_output.cc
 * @brief   
 */

#include "input/type_output.hh"
#include "input/type_repository.hh"
#include "input/type_generic.hh"
#include "input/type_tuple.hh"
#include "system/system.hh"
#include <boost/algorithm/string/replace.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/concepts.hpp>
#include <boost/iostreams/operations.hpp> // put


#include <string>
#include <limits>
#include <boost/regex.hpp>

namespace Input {
namespace Type {

using namespace std;


/*******************************************************************
 * implementation of OutputBase
 */

OutputBase::~OutputBase() {}



OutputBase::OutputBase()
{
    TypeBase::lazy_finish();
}



void OutputBase::get_integer_bounds(Integer integer, int &lower , int &upper ) {
    lower = integer.lower_bound_;
    upper = integer.upper_bound_;
}



void OutputBase::get_double_bounds(Double dbl, double &lower , double &upper ) {
    lower = dbl.lower_bound_;
    upper = dbl.upper_bound_;
}



void OutputBase::get_array_sizes(Array array, unsigned int &lower , unsigned int &upper ) {
	lower = array.data_->lower_bound_;
	upper = array.data_->upper_bound_;
}



void OutputBase::get_array_type(Array array, std::shared_ptr<TypeBase> &arr_type) {
    arr_type = array.data_->type_of_values_;
}



const string & OutputBase::get_record_description(const Record *rec) {
    return rec->data_->description_;
}



const string & OutputBase::get_abstract_description(const Abstract *a_rec) {
    return a_rec->child_data_->description_;
}



void OutputBase::get_parent_vec(Record rec, std::vector< std::shared_ptr<Abstract> > &parent_vec) {
	parent_vec = rec.data_->parent_vec_;
}



void OutputBase::get_default(Record::KeyIter it, string &type, string &value) {
	value = it->default_.value_;
	if ( it->default_.is_obligatory() ) {
		type = "obligatory";
	} else if ( it->default_.is_optional() ) {
		type = "optional";
	} else if ( it->default_.has_value_at_read_time() ) {
		type = "value at read time";
	} else {
		type = "value at declaration";
	}

}


const string & OutputBase::get_selection_description(const Selection *sel) {
    return sel->data_->description_;
}


void OutputBase::get_adhoc_parent_name(const AdHocAbstract *a_rec, string &parent_name) {
	parent_name = a_rec->ancestor_.type_name();
}


Abstract::ChildDataIter OutputBase::get_adhoc_parent_data(const AdHocAbstract *a_rec) {
	return a_rec->ancestor_.child_data_->list_of_childs.begin();
}




void OutputBase::print_base(ostream& stream, const TypeBase *type) {

	if (typeid(*type) == typeid(Type::Record) || typeid(*type) == typeid(Type::Tuple)) {
		print_impl(stream, static_cast<const Type::Record *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Array)) {
		print_impl(stream, static_cast<const Type::Array *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Abstract)) {
			print_impl(stream, static_cast<const Type::Abstract *>(type) );
	} else
	if (typeid(*type) == typeid(Type::AdHocAbstract)) {
		print_impl(stream, static_cast<const Type::AdHocAbstract *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Selection)) {
		print_impl(stream, static_cast<const Type::Selection *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Integer)) {
		print_impl(stream, static_cast<const Type::Integer *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Double)) {
		print_impl(stream, static_cast<const Type::Double *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Bool)) {
		print_impl(stream, static_cast<const Type::Bool *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Parameter)) {
		print_impl(stream, static_cast<const Type::Parameter *>(type) );
	} else {
		const Type::FileName * file_name_type = dynamic_cast<const Type::FileName *>(type);
        if (file_name_type != NULL ) {
        	print_impl(stream, file_name_type );
        	return;
        }

		const Type::String * string_type = dynamic_cast<const Type::String *>(type);
        if (string_type != NULL ) {
        	print_impl(stream, string_type );
        	return;
        }

        // default -> error
        THROW( Type::ExcUnknownDescendant() << Type::EI_TypeName(typeid(type).name()) );
	}
}


void OutputBase::write_default_value(std::ostream& stream, Default dft) {
	if (dft.is_obligatory() || dft.is_optional()) {
		stream << "<" << dft.value() << ">";
	} else {
		stream << "\"" << dft.value() << "\"";
	}
}

void OutputBase::write_description(std::ostream& stream, const string& str,
        unsigned int padding, unsigned int hash_count) {
	string s = str;
	boost::replace_all(s, "\\$", "$");

    boost::tokenizer<boost::char_separator<char> > line_tokenizer(s, boost::char_separator<char>("\n"));
    boost::tokenizer<boost::char_separator<char> >::iterator it;

    // For every \n add padding at beginning of the next line.
    for(it = line_tokenizer.begin(); it != line_tokenizer.end(); ++it) {
        stream << endl;
        stream << setw(padding) << "";
        stream << std::setfill('#') << setw(hash_count) << "" << std::setfill(' ') << " " << *it;
    }
}


void OutputBase::clear_processed_types() {
	processed_types_hash_.clear();
	full_hash_ = 0;
}


bool OutputBase::was_written(std::size_t hash)
{
    bool in_set = ( processed_types_hash_.find(hash) != processed_types_hash_.end() );
    if (! in_set) processed_types_hash_.insert(hash);
    return in_set;
}

void OutputBase::get_attr_and_param_data(const TypeBase *type, TypeBase::attribute_map &attr_map,
    		TypeBase::TypeHash &generic_type_hash, TypeBase::json_string &parameter_map_to_json) {
	attr_map = *( type->attributes_.get() );
	generic_type_hash = type->generic_type_hash_;
	if (type->parameter_map_.size()) {
		parameter_map_to_json = type->print_parameter_map_to_json( type->parameter_map_ );
	}
}




/*******************************************************************
 * implementation of OutputText
 */

ostream& OutputText::print(ostream& stream) {
	doc_type_ = full_record;
	clear_processed_types();

	print_base(stream, type_);
	return stream;
}

void OutputText::print_impl(ostream& stream, const Record *type) {
	ASSERT(type->is_finished())(type->type_name())(type->class_name()).warning("Printing documentation of unfinished type.");
	switch (doc_type_) {
	case key_record:
		stream << "" << type->class_name() << " '" << type->type_name() << "' (" << type->size() << " keys).";
		break;
	case full_record:
		TypeBase::TypeHash hash=type->content_hash();
		if (! was_written(hash)) {
			// header
			stream << endl;
			stream << "" << type->class_name() << " '" << type->type_name() << "'";
			// parent record
			/*std::shared_ptr<Abstract> parent_ptr;
			get_parent_ptr(*type, parent_ptr);
			if (parent_ptr) {
				stream << ", implementation of " << parent_ptr->type_name();
			}*/

			// reducible to key
			Record::KeyIter key_it = type->auto_conversion_key_iter();
			if (key_it != type->end()) {
				stream << ", reducible to key '" << key_it->key_ << "'";
			}
			stream << "" << " (" << type->size() << " keys).";
		    write_description(stream, OutputBase::get_record_description(type), 0);
		    stream << endl;
		    stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
		    // keys
		    doc_type_ = key_record;
		    for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
		    	size_setw_ = it->key_.size() + 3;
		        stream << setw(padding_size) << "" << it->key_ << " = ";
		        write_default_value(stream, it->default_);
		        stream << endl;
		        stream << setw(padding_size + size_setw_) << "" <<"#### is ";
		        print_base(stream, it->type_.get());
		        write_description(stream, it->description_, padding_size+size_setw_);
		        stream << endl;
		    }
		    stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type->type_name() << endl;
		    // Full documentation of embedded record types.
		    doc_type_ = full_record;
		}
		break;
	}
}
void OutputText::print_impl(ostream& stream, const Array *type) {
	std::shared_ptr<TypeBase> array_type;
	get_array_type(*type, array_type);
	switch (doc_type_) {
	case key_record:
		unsigned int lower_size, upper_size;
		get_array_sizes(*type, lower_size, upper_size);
		stream << "Array, size limits: [" << lower_size << ", " << upper_size << "] of type: " << endl;
		stream << setw(padding_size + size_setw_) << "" << "#### ";
		print_base(stream, array_type.get());
		break;
	case full_record:
		print_base(stream, array_type.get());
		break;
	}
}


void OutputText::print_impl(ostream& stream, const Abstract *type) {
	// Print documentation of abstract record
	switch (doc_type_) {
	case key_record:
		stream << "Abstract '" << type->type_name() << "' with "<< type->child_size() << " descendants.";
		break;
	case full_record:
		TypeBase::TypeHash hash=type->content_hash();
		if (! was_written(hash) ) {
            // header
            stream << endl;
            stream << "" << "Abstract '" << type->type_name() << "' with " << type->child_size() << " descendants.";
            write_description(stream, OutputBase::get_abstract_description( type ), 0);
            stream << endl;
            stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
            // descendants
            doc_type_ = key_record;
            for (Abstract::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
            	size_setw_ = 0;
                stream << setw(padding_size) << "";
                stream << "" << "Record '" << (*it).type_name() << "'";
                write_description(stream, OutputBase::get_record_description( &(*it) ), padding_size+size_setw_);
                stream << endl;
            }
            stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type->type_name() << endl;
            // Full documentation of embedded record types.
            doc_type_ = full_record;
        }
		break;
	}
}


void OutputText::print_impl(ostream& stream, const AdHocAbstract *type) {
	// Print documentation of adhoc abstract record
	if (doc_type_ == key_record) {
		stream << "AdHocAbstract";
		/*string parent_name;
		get_adhoc_parent_name(type, parent_name);
		stream << "AdHocAbstract" << endl;
		stream << setw(padding_size + size_setw_) << "";
		stream << "#### Derived from Abstract '" << parent_name << "', ";
		stream << "added Records: ";
		{
			Abstract::ChildDataIter parent_it = get_adhoc_parent_data(type);
			bool add_comma = false;
			for (Abstract::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
				if ((*it).type_name() == (*parent_it).type_name()) {
					++parent_it;
				} else {
					if (add_comma) stream << ", ";
					else add_comma = true;
					stream << "'" << (*it).type_name() << "'";
				}
			}
		}*/
	}
}
void OutputText::print_impl(ostream& stream, const Selection *type) {
	ASSERT(type->is_finished())(type->type_name()).warning("Printing documentation of unfinished Input::Type::Selection.");
	switch (doc_type_) {
	case key_record:
		stream << "Selection '" << type->type_name() << "' of " << type->size() << " values.";
		break;
	case full_record:
		TypeBase::TypeHash hash=type->content_hash();
		if (! was_written(hash) ) {
			stream << endl << "Selection '" << type->type_name() << "' of " << type->size() << " values.";
			write_description(stream, OutputBase::get_selection_description( type ), 0);
			stream << endl;
		    stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
		    // keys
		    for (Selection::keys_const_iterator it = type->begin(); it != type->end(); ++it) {
		        stream << setw(padding_size) << "" << it->key_ << " = " << it->value;
		        if (it->description_ != "") {
		        	stream << endl;
		        	stream << setw(padding_size + it->key_.size() + 3) << "" << "# " << it->description_ << "";
		        }
		        stream << endl;
		    }
		    stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type->type_name() << endl;
		}
		break;
	}
}
void OutputText::print_impl(ostream& stream, const Integer *type) {
	if (doc_type_ == key_record) {
		int lower_bound, upper_bound;
		get_integer_bounds(*type, lower_bound, upper_bound);
		stream << "Integer in [" << lower_bound << ", " << upper_bound << "]";
	}
}
void OutputText::print_impl(ostream& stream, const Double *type) {
	if (doc_type_ == key_record) {
		double lower_bound, upper_bound;
		get_double_bounds(*type, lower_bound, upper_bound);
		stream << "Double in [" << lower_bound << ", " << upper_bound << "]";
	}
}
void OutputText::print_impl(ostream& stream, const Bool *type) {
	if (doc_type_ == key_record) {
		stream << "Bool";
	}
}
void OutputText::print_impl(ostream& stream, const String *type) {
	if (doc_type_ == key_record) {
		stream << "String (generic)";
	}
}
void OutputText::print_impl(ostream& stream, const FileName *type) {
	if (doc_type_ == key_record) {
		stream << "FileName of ";
		switch (type->get_file_type()) {
		case ::FilePath::input_file:
			stream << "input file";
			break;
		case ::FilePath::output_file:
			stream << "output file";
			break;
		default:
			stream << "file with unknown type";
			break;
		}
	}
}



void OutputText::print_impl(ostream& stream, const Parameter *type) {
	ASSERT_DBG(false).error("Parameter appears in the IST. Check where Instance is missing.");
}




/*******************************************************************
 * implementation of OutputJSONMachine
 */

OutputJSONMachine::OutputJSONMachine(RevNumData rev_num_data) 
: OutputBase()
{    
    rev_num_data_ = rev_num_data;

    format_head="{ \"version\" :";
    format_inner=",\n\"ist_nodes\" : [\n";
    format_full_hash="{}],\n";
    format_tail="}\n";   
}



OutputJSONMachine::OutputJSONMachine(const Record &root_type, RevNumData rev_num_data) 
: OutputJSONMachine(rev_num_data)
{
    root_type_ =  root_type;
}


std::string OutputJSONMachine::escape_description(std::string desc) {
	static std::vector< std::pair<boost::regex, std::string> > rewrite_rules = {
	        // replace single slash with two slashes
			{boost::regex("\\\\"), "\\\\\\\\"},
	        // replace quote with slash quote
			{boost::regex("\\\""), "\\\\\""},
	        // replace special chars with escaped slash + special chars
			{boost::regex("\\n"), "\\\\n"},
			{boost::regex("\\t"), "\\\\t"},
			{boost::regex("\\r"), "\\\\r"}
	};


    std::string tmp = std::string(desc);

    for (auto rewrite_rule : rewrite_rules) {
        tmp = boost::regex_replace(tmp, rewrite_rule.first, rewrite_rule.second);
    }

    return tmp;
}


ostream& OutputJSONMachine::print(ostream& stream) {
	doc_type_ = full_record;
	clear_processed_types();

	stream << format_head;
	print_program_info(stream);
	stream << format_inner;

    print_base( stream, &root_type_);
	for (Input::TypeRepository<Selection>::TypeRepositoryMapIter it = Input::TypeRepository<Selection>::get_instance().begin();
			it != Input::TypeRepository<Selection>::get_instance().end(); ++it) {
		print_base( stream, it->second.get() );
	}
	for (Input::TypeRepository<Abstract>::TypeRepositoryMapIter it = Input::TypeRepository<Abstract>::get_instance().begin();
			it != Input::TypeRepository<Abstract>::get_instance().end(); ++it) {
		print_base( stream, it->second.get() );
	}
	for (Input::TypeRepository<Record>::TypeRepositoryMapIter it = Input::TypeRepository<Record>::get_instance().begin();
			it != Input::TypeRepository<Record>::get_instance().end(); ++it) {
		print_base( stream, it->second.get() );
	}

	stream << format_full_hash;
	print_full_hash(stream);
	stream << format_tail;
	return stream;
}

std::string print_attributes(TypeBase::attribute_map attribute_map) {
    stringstream stream;

    if (attribute_map.size() == 0) return "\"attributes\" : {}";

    stream << "\"attributes\" : {" << endl; // print map of attributes
    for (auto it=attribute_map.begin(); it!=attribute_map.end(); ++it) {
        if (it != attribute_map.begin()) {
            stream << "," << endl;
        }
        stream << "\"" << it->first << "\" : " << it->second;
    }
    stream << endl << "}";

    return stream.str();
}


void OutputJSONMachine::print_type_header(ostream &stream, const TypeBase *type) {
    stream << "{" << endl;
    stream << "\"id\" : " << type->hash_str() << "," << endl;
    stream << "\"input_type\" : \"" + type->class_name() + "\"," << endl;
    stream << "\"name\" : \"" << type->type_name() << "\"," << endl;

    // write data of parameters and attributes
    TypeBase::attribute_map attr_map;
    TypeBase::TypeHash generic_type_hash;
    TypeBase::json_string parameter_map_to_json;
    get_attr_and_param_data(type, attr_map, generic_type_hash, parameter_map_to_json);
	// print hash of generic type and parameters into separate keys
	if (generic_type_hash) { // print hash of generic type into separate keys
		stream << "\"generic_type\" : " << TypeBase::hash_str(generic_type_hash) << "," << endl;
	}
	if (parameter_map_to_json.size()) { // parameters into separate keys
		stream << "\"parameters\" : " << parameter_map_to_json << "," << endl;
	}
	stream << print_attributes(attr_map);
}


void OutputJSONMachine::print_impl(ostream& stream, const Record *type) {

	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    print_type_header(stream, type);
    stream << "," << endl << "\"description\" : \"" <<
            escape_description( OutputBase::get_record_description(type) ) << "\"," << endl;

    // parent records, implemented abstracts
    std::vector< std::shared_ptr<Abstract> > parent_vec;
    get_parent_vec(*type, parent_vec);
    if (parent_vec.size()) {
        stream << "\"implements\" : [ ";
        bool add_comma = false;
        for (auto &parent : parent_vec) {
        	if (add_comma) stream << ", ";
        	else add_comma = true;
            stream << parent->hash_str();
        }
        stream << " ]," << endl;
    }

    // reducible to key
    Record::KeyIter key_it = type->auto_conversion_key_iter();
    if (key_it != type->end()) {
        stream << "\"reducible_to_key\" : \"" << key_it->key_ << "\"," << endl;
    }

    stream << "\"keys\" : [" << endl;

    doc_type_ = key_record;
    for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
        string dft_type, dft_value;
        get_default(it, dft_type, dft_value);
        if (dft_type != "value at declaration") 
            dft_value = "\"" + escape_description(dft_value) + "\"";

        if (it != type->begin()) {
            stream << "," << endl;
        }
        stream << "{ \"key\" : \"" << it->key_ << "\"," << endl;
        stream << "\"description\" : \"" <<
                escape_description(it->description_) << "\"," << endl;
        stream << "\"default\" : { "
                <<"\"type\" : \"" << dft_type << "\"," << endl
                <<"\"value\" : " << dft_value << " }," << endl;
        stream << "\"type\" : " << it->type_->hash_str() << "," << endl;
        stream << print_attributes(it->attributes);
        stream << "}";
    }

    stream << "]" << endl;
    stream << "},";

    // Full documentation of embedded record types.
    doc_type_ = full_record;
	for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
		print_base(stream, it->type_.get());
	}

    boost::hash_combine(full_hash_, hash);
}



void OutputJSONMachine::print_impl(ostream& stream, const Array *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    unsigned int lower_size, upper_size;
	std::shared_ptr<TypeBase> array_type;

    get_array_sizes(*type, lower_size, upper_size);
	get_array_type(*type, array_type);

	print_type_header(stream, type);
	stream << "," << endl;
	stream << "\"range\" : [" << lower_size << ", " << upper_size << "]," << endl;
	stream << "\"subtype\" : " << array_type->hash_str()  << endl;
	stream << "}," << endl;

	print_base(stream, array_type.get());

	boost::hash_combine(full_hash_, hash);
}



void OutputJSONMachine::print_impl(ostream& stream, const Abstract *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    print_type_header(stream, type);
    stream << "," << endl;
    stream << "\"description\" : \"" <<
            escape_description( OutputBase::get_abstract_description(type)) << "\"," << endl;

    print_abstract_record_keys(stream, type);
    stream << "},";

    for (Abstract::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
    	print_base(stream, &*it);
    }

    boost::hash_combine(full_hash_, hash);
}


void OutputJSONMachine::print_impl(ostream& stream, const AdHocAbstract *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

	string parent_name;
	get_adhoc_parent_name(type, parent_name);
	print_type_header(stream, type);
	stream << "," << endl;
    print_abstract_record_keys(stream, dynamic_cast<const Type::Abstract *>(type));
    stream << "},";

    for (Abstract::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
    	print_base(stream, &*it);
    }

    boost::hash_combine(full_hash_, hash);
}



void OutputJSONMachine::print_abstract_record_keys(ostream& stream, const Abstract *type) {

    // Print documentation of abstract record
    const Record * desc = type->get_default_descendant();

    // default descendant
    if (desc) {
        stream << "\"default_descendant\" : " << desc->hash_str()  << "," << endl;
    }
    stream << "\"implementations\" : [" << endl;
    for (Abstract::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
        if (it != type->begin_child_data())
            stream << "," << endl;
        stream << "" << it->hash_str() << "";
    }
    stream << "]";

}





void OutputJSONMachine::print_impl(ostream& stream, const Selection *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    print_type_header(stream, type);
    stream << "," << endl;
	stream << "\"description\" : \"" <<
	        escape_description(OutputBase::get_selection_description(type)) << "\"," << endl;

	stream << "\"values\" : [" << endl;
	for (Selection::keys_const_iterator it = type->begin(); it != type->end(); ++it) {
		if (it != type->begin()) {
			stream << "," << endl;
		}
		stream << "{ \"name\" : \"" << it->key_ << "\"," << endl;
		stream << "\"description\" : \"" << escape_description(it->description_) << "\"," << endl;
		stream << print_attributes(it->attributes_) << "," << endl;
		stream << "}";
	}

	stream << "]" << endl;
	stream << "},";

	boost::hash_combine(full_hash_, hash);
}


void OutputJSONMachine::print_impl(ostream& stream, const Integer *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    int lower, upper;
	get_integer_bounds(*type, lower, upper);

	print_type_header(stream, type);
	stream << "," << endl;
	stream << "\"range\" : [" << lower << ", " << upper << "]" << endl;
	stream << "},";

	boost::hash_combine(full_hash_, hash);
}


void OutputJSONMachine::print_impl(ostream& stream, const Double *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    double lower, upper;
	get_double_bounds(*type, lower, upper);

	print_type_header(stream, type);
	stream << "," << endl;
    stream << "\"range\" : [" << lower << ", " << upper << "]" << endl;
	stream << "},";

	boost::hash_combine(full_hash_, hash);
}


void OutputJSONMachine::print_impl(ostream& stream, const Bool *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    print_type_header(stream, type);
	stream << "},";

	boost::hash_combine(full_hash_, hash);
}


void OutputJSONMachine::print_impl(ostream& stream, const String *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    print_type_header(stream, type);
	stream << "},";

	boost::hash_combine(full_hash_, hash);
}


void OutputJSONMachine::print_impl(ostream& stream, const FileName *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    print_type_header(stream, type);
    stream << "," << endl;
	stream << "\"input_type\" : \"FileName\"," << endl;
	stream << "\"file_mode\" : \"";
	switch (type->get_file_type()) {
	case ::FilePath::input_file:
		stream << "input\"";
		break;
	case ::FilePath::output_file:
		stream << "output\"";
		break;
	}

	stream << endl << "},";

	boost::hash_combine(full_hash_, hash);
}



void OutputJSONMachine::print_impl(ostream& stream, const Parameter *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    print_type_header(stream, type);
    stream << endl << "},";

	boost::hash_combine(full_hash_, hash);
}



void OutputJSONMachine::print_program_info(ostream& stream) {
	string build_date = string(__DATE__) + ", " + string(__TIME__);

    stream << "{" << endl;
    stream << "\"flow123d_commit\" : \"" << rev_num_data_.revision << "\"," << endl;
    stream << "\"flow123d_version\" : \"" << rev_num_data_.version << "\"," << endl;
    stream << "\"date\" : \"" << build_date << "\"" << endl;
	stream << "}";
}



void OutputJSONMachine::print_full_hash(ostream& stream) {
	stream << "\"IST_hash\" : " << TypeBase::hash_str( full_hash_ ) << endl;
}






std::ostream& operator<<(std::ostream& stream, OutputText type_output) {
    return type_output.print(stream) << endl;
}


std::ostream& operator<<(std::ostream& stream, OutputJSONMachine type_output) {
    return type_output.print(stream) << endl;
}


} // closing namespace Type
} // closing namespace Input
