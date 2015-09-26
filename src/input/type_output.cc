/**
 * type_output.cc
 */

#include "input/type_output.hh"
#include "input/type_repository.hh"
#include "rev_num.h"
#include <boost/algorithm/string/replace.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/concepts.hpp>
#include <boost/iostreams/operations.hpp> // put


#include <string>
#include <limits>

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



void OutputBase::get_array_type(Array array, boost::shared_ptr<TypeBase> &arr_type) {
    arr_type = array.data_->type_of_values_;
}



const string & OutputBase::get_record_description(const Record *rec) {
    return rec->data_->description_;
}



const string & OutputBase::get_abstract_description(const AbstractRecord *a_rec) {
    return a_rec->child_data_->description_;
}



void OutputBase::get_parent_ptr(Record rec, boost::shared_ptr<AbstractRecord> &parent_ptr) {
	if (rec.data_->parent_vec_.size()) {
		// temporary solution, if we need this method in new input, it must return vector of parents
		parent_ptr = rec.data_->parent_vec_[0];
	} else {
		parent_ptr = boost::shared_ptr<AbstractRecord>();
	}
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


AbstractRecord::ChildDataIter OutputBase::get_adhoc_parent_data(const AdHocAbstractRecord *a_rec) {
	return a_rec->parent_data_->list_of_childs.begin();
}



const string & OutputBase::get_adhoc_parent_name(const AdHocAbstractRecord *a_rec) {
	return a_rec->parent_name_;
}



void OutputBase::print_base(ostream& stream, const TypeBase *type) {

	if (typeid(*type) == typeid(Type::Record)) {
		print_impl(stream, static_cast<const Type::Record *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Array)) {
		print_impl(stream, static_cast<const Type::Array *>(type) );
	} else
	if (typeid(*type) == typeid(Type::AbstractRecord)) {
			print_impl(stream, static_cast<const Type::AbstractRecord *>(type) );
	} else
	if (typeid(*type) == typeid(Type::AdHocAbstractRecord)) {
		print_impl(stream, static_cast<const Type::AdHocAbstractRecord *>(type) );
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
        xprintf(Err,"Unknown descendant of TypeBase class, name: %s\n", typeid(type).name());
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
	if (! type->is_finished()) {
		xprintf(Warn, "Printing documentation of unfinished Input::Type::Record!\n");
	}
	switch (doc_type_) {
	case key_record:
		stream << "" << "Record '" << type->type_name() << "' (" << type->size() << " keys).";
		break;
	case full_record:
		TypeBase::TypeHash hash=type->content_hash();
		if (! was_written(hash)) {
			// header
			stream << endl;
			stream << "" << "Record '" << type->type_name() << "'";
			// parent record
			boost::shared_ptr<AbstractRecord> parent_ptr;
			get_parent_ptr(*type, parent_ptr);
			if (parent_ptr) {
				stream << ", implementation of " << parent_ptr->type_name();
			}
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
	boost::shared_ptr<TypeBase> array_type;
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
void OutputText::print_impl(ostream& stream, const AbstractRecord *type) {
	// Print documentation of abstract record
	switch (doc_type_) {
	case key_record:
		stream << "AbstractRecord '" << type->type_name() << "' with "<< type->child_size() << " descendants.";
		break;
	case full_record:
		TypeBase::TypeHash hash=type->content_hash();
		if (! was_written(hash) ) {
            // header
            stream << endl;
            stream << "" << "AbstractRecord '" << type->type_name() << "' with " << type->child_size() << " descendants.";
            write_description(stream, OutputBase::get_abstract_description( type ), 0);
            stream << endl;
            stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
            // descendants
            doc_type_ = key_record;
            for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
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
void OutputText::print_impl(ostream& stream, const AdHocAbstractRecord *type) {
	// Print documentation of adhoc abstract record
	if (doc_type_ == key_record) {
		stream << "AdHocAbstractRecord" << endl;
		stream << setw(padding_size + size_setw_) << "";
		stream << "#### Derived from AbstractRecord '" << get_adhoc_parent_name(type) << "', ";
		stream << "added Records: ";
		{
			AbstractRecord::ChildDataIter parent_it = get_adhoc_parent_data(type);
			bool add_comma = false;
			for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
				if ((*it).type_name() == (*parent_it).type_name()) {
					++parent_it;
				} else {
					if (add_comma) stream << ", ";
					else add_comma = true;
					stream << "'" << (*it).type_name() << "'";
				}
			}
		}
	}
}
void OutputText::print_impl(ostream& stream, const Selection *type) {
	if (! type->is_finished()) {
		xprintf(Warn, "Printing documentation of unfinished Input::Type::Selection!\n");
	}
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






/*******************************************************************
 * implementation of OutputJSONMachine
 */


std::string OutputJSONMachine::format_hash( TypeBase::TypeHash hash) {
    stringstream ss;
    ss << std::hex << hash;
    return ss.str();
}


std::string OutputJSONMachine::escape_description(std::string desc) {
    static OutputJSONMachine::RewriteRule rewrite_rules[] = {
        // replace single slash with two slashes 
        OutputJSONMachine::RewriteRule (boost::regex("\\\\"), "\\\\\\\\"),
        // replace quote with slash quote
        OutputJSONMachine::RewriteRule (boost::regex("\\\""), "\\\\\""),
        // replace special chars with escaped slash + special chars
        OutputJSONMachine::RewriteRule (boost::regex("\\n"), "\\\\n"),
        OutputJSONMachine::RewriteRule (boost::regex("\\t"), "\\\\t"),
        OutputJSONMachine::RewriteRule (boost::regex("\\r"), "\\\\r")
    };


    std::string tmp = std::string(desc);

    for (OutputJSONMachine::RewriteRule rewrite_rule : rewrite_rules) {
        tmp = boost::regex_replace(tmp, rewrite_rule.search, rewrite_rule.replacement);
    }

    return tmp;
}


ostream& OutputJSONMachine::print(ostream& stream) {
	doc_type_ = full_record;
	clear_processed_types();

	stream << format_head;
	print_program_info(stream);
	stream << format_inner;

	for (Input::TypeRepository<Selection>::TypeRepositoryMapIter it = Input::TypeRepository<Selection>::get_instance().begin();
			it != Input::TypeRepository<Selection>::get_instance().end(); ++it) {
		print_base( stream, it->second.get() );
	}
	for (Input::TypeRepository<AbstractRecord>::TypeRepositoryMapIter it = Input::TypeRepository<AbstractRecord>::get_instance().begin();
			it != Input::TypeRepository<AbstractRecord>::get_instance().end(); ++it) {
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


void OutputJSONMachine::print_impl(ostream& stream, const Record *type) {

	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"Record\"," << endl;
    stream << "\"type_name\" : \"" << type->type_name() << "\"," << endl;
    type->write_attributes(stream);
    stream << "," << endl << endl;
    stream << "\"description\" : \"" <<
            escape_description( OutputBase::get_record_description(type) ) << "\"," << endl;

    // parent records, implemented abstracts
    boost::shared_ptr<AbstractRecord> parent_ptr;
    get_parent_ptr(*type, parent_ptr);
    if (parent_ptr) {
        stream << "\"implements\" : [ \"" << format_hash(parent_ptr->content_hash()) << "\" ]," << endl;
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

        if (it != type->begin()) {
            stream << "," << endl;
        }
        stream << "{ \"key\" : \"" << it->key_ << "\"," << endl;
        stream << "\"description\" : \"" <<
                escape_description(it->description_) << "\"," << endl;
        stream << "\"default\" : { "
                <<"\"type\" : \"" << dft_type << "\"," << endl
                <<"\"value\" : \"" << escape_description(dft_value) << "\" }," << endl;
        stream << "\"type\" : \"" << format_hash(it->type_->content_hash()) << "\"" << endl;
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
	boost::shared_ptr<TypeBase> array_type;

    get_array_sizes(*type, lower_size, upper_size);
	get_array_type(*type, array_type);

	stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"Array\"," << endl;
	stream << "\"range\" : [" << lower_size << ", " << upper_size << "]," << endl;
	stream << "\"subtype\" : \"" << format_hash(array_type->content_hash()) << "\"," << endl;
	type->write_attributes(stream);
	stream << endl;
	stream << "}," << endl;

	print_base(stream, array_type.get());

	boost::hash_combine(full_hash_, hash);
}



void OutputJSONMachine::print_impl(ostream& stream, const AbstractRecord *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"AbstractRecord\"," << endl;
    stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
    type->write_attributes(stream);
    stream << "," << endl;
    stream << "\"description\" : \"" <<
            escape_description( OutputBase::get_abstract_description(type)) << "\"," << endl;

    print_abstract_record_keys(stream, type);
    stream << "},";

    for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
    	print_base(stream, &*it);
    }

    boost::hash_combine(full_hash_, hash);
}


void OutputJSONMachine::print_impl(ostream& stream, const AdHocAbstractRecord *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"AdHocAbstractRecord\"," << endl;
    stream << "\"parent\" : \"" << get_adhoc_parent_name(type) << "\"," << endl;
    type->write_attributes(stream);
    stream << "," << endl;

    print_abstract_record_keys(stream, dynamic_cast<const Type::AbstractRecord *>(type));
    stream << "},";

    for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
    	print_base(stream, &*it);
    }

    boost::hash_combine(full_hash_, hash);
}



void OutputJSONMachine::print_abstract_record_keys(ostream& stream, const AbstractRecord *type) {

    // Print documentation of abstract record
    const Record * desc = type->get_default_descendant();

    // default descendant
    if (desc) {
        stream << "\"default_descendant\" : \"" << format_hash(desc->content_hash())  << "\"," << endl;
    }
    stream << "\"implementations\" : [" << endl;
    for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
        if (it != type->begin_child_data()) {
            stream << ",\n" << endl;
        }

        stream << "\"" << format_hash(it->content_hash()) << "\"";
    }
    stream << "]";

}





void OutputJSONMachine::print_impl(ostream& stream, const Selection *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

	stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"Selection\"," << endl;
    stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	type->write_attributes(stream);
    stream << "," << endl;
	stream << "\"description\" : \"" <<
	        escape_description(OutputBase::get_selection_description(type)) << "\"," << endl;

	stream << "\"values\" : [" << endl;

	for (Selection::keys_const_iterator it = type->begin(); it != type->end(); ++it) {
		if (it != type->begin()) {
			stream << "," << endl;
		}
		stream << "{ \"name\" : \"" << it->key_ << "\"," << endl
		       << "\"description\" : \"" << escape_description(it->description_) << "\" }";
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

	stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"Integer\"," << endl;
    stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	type->write_attributes(stream);
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

	stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"Double\"," << endl;
    stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	type->write_attributes(stream);
    stream << "," << endl;

    stream << "\"range\" : [" << lower << ", " << upper << "]" << endl;
	stream << "},";

	boost::hash_combine(full_hash_, hash);
}


void OutputJSONMachine::print_impl(ostream& stream, const Bool *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"Bool\"," << endl;
    stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	type->write_attributes(stream);
    stream << endl;
	stream << "},";

	boost::hash_combine(full_hash_, hash);
}


void OutputJSONMachine::print_impl(ostream& stream, const String *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"String\"," << endl;
    stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	type->write_attributes(stream);
    stream << endl;
	stream << "},";

	boost::hash_combine(full_hash_, hash);
}


void OutputJSONMachine::print_impl(ostream& stream, const FileName *type) {
	TypeBase::TypeHash hash=type->content_hash();
    if (was_written(hash)) return;

    stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
	stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	type->write_attributes(stream);
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



void OutputJSONMachine::print_program_info(ostream& stream) {
    string revision(_GIT_REVISION_);
	string version(_VERSION_NAME_);
	string build_date = string(__DATE__) + ", " + string(__TIME__);

    stream << "{" << endl;
    stream << "\"flow123d_commit\" : \"" << revision << "\"," << endl;
    stream << "\"flow123d_version\" : \"" << version << "\"," << endl;
    stream << "\"date\" : \"" << build_date << "\"" << endl;
	stream << "}";
}



void OutputJSONMachine::print_full_hash(ostream& stream) {
	stream << "\"IST_hash\" : \"" << format_hash(full_hash_) << "\"" << endl;
}






std::ostream& operator<<(std::ostream& stream, OutputText type_output) {
    return type_output.print(stream) << endl;
}


std::ostream& operator<<(std::ostream& stream, OutputJSONMachine type_output) {
    return type_output.print(stream) << endl;
}


} // closing namespace Type
} // closing namespace Input
