/**
 * type_output.cc
 */

#include "input/type_output.hh"
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



OutputBase::OutputBase(const TypeBase *type, unsigned int depth)
: type_(type), depth_(depth)
{
    TypeBase::lazy_finish();
}



ostream&  OutputBase::print(ostream& stream) {
	doc_type_ = full_record;
	doc_flags_.clear();

	stream << format_head;
	print(stream, type_, depth_);
	stream << format_tail;
	return stream;
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



void OutputBase::get_parent_vec(Record rec, std::vector< boost::shared_ptr<AbstractRecord> > &parent_vec) {
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


const string & OutputBase::get_adhoc_parent_name(const AdHocAbstractRecord *a_rec) {
	return a_rec->parent_name_;
}



void OutputBase::print(ostream& stream, const TypeBase *type, unsigned int depth) {
	if (typeid(*type) == typeid(Type::Record)) {
		print_impl(stream, static_cast<const Type::Record *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::Array)) {
		print_impl(stream, static_cast<const Type::Array *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::AbstractRecord)) {
			print_impl(stream, static_cast<const Type::AbstractRecord *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::AdHocAbstractRecord)) {
		print_impl(stream, static_cast<const Type::AdHocAbstractRecord *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::Selection)) {
		print_impl(stream, static_cast<const Type::Selection *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::Integer)) {
		print_impl(stream, static_cast<const Type::Integer *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::Double)) {
		print_impl(stream, static_cast<const Type::Double *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::Bool)) {
		print_impl(stream, static_cast<const Type::Bool *>(type), depth );
	} else {
		const Type::FileName * file_name_type = dynamic_cast<const Type::FileName *>(type);
        if (file_name_type != NULL ) {
        	print_impl(stream, file_name_type, depth );
        	return;
        }

		const Type::String * string_type = dynamic_cast<const Type::String *>(type);
        if (string_type != NULL ) {
        	print_impl(stream, string_type, depth );
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




/*******************************************************************
 * implementation of OutputBase::ProcessedTypes
 */


OutputBase::ProcessedTypes::~ProcessedTypes() {}


void OutputBase::ProcessedTypes::clear() {
	output_hash.clear();
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


void OutputJSONMachine::print_impl(ostream& stream, const Record *type, unsigned int depth) {

	TypeBase::TypeHash hash=type->content_hash();
    if (doc_flags_.was_written(hash)) return;

    stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"Record\"," << endl;
    stream << "\"type_name\" : \"" << type->type_name() << "\"," << endl;
    type->write_attributes(stream);
    stream << "," << endl << endl;
    stream << "\"description\" : \"" <<
            escape_description( OutputBase::get_record_description(type) ) << "\"," << endl;

    // parent records, implemented abstracts
    std::vector< boost::shared_ptr<AbstractRecord> > parent_vec;
    get_parent_vec(*type, parent_vec);
    if (parent_vec.size()) {
        stream << "\"implements\" : [ ";
        bool add_comma = false;
        for (auto &parent : parent_vec) {
        	if (add_comma) stream << ", ";
        	else add_comma = true;
            stream << "\"" << format_hash(parent->content_hash()) << "\"";
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
    if (depth_ == 0 || depth_ > depth) {
        for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
            print(stream, it->type_.get(), depth+1);
        }
    }
}



void OutputJSONMachine::print_impl(ostream& stream, const Array *type, unsigned int depth) {
	TypeBase::TypeHash hash=type->content_hash();
    if (doc_flags_.was_written(hash)) return;

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

	print(stream, array_type.get() ,depth+1);
}



void OutputJSONMachine::print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth) {
	TypeBase::TypeHash hash=type->content_hash();
    if (doc_flags_.was_written(hash)) return;

    stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"AbstractRecord\"," << endl;
    stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
    type->write_attributes(stream);
    stream << "," << endl;
    stream << "\"description\" : \"" <<
            escape_description( OutputBase::get_abstract_description(type)) << "\"," << endl;

    print_abstract_record_keys(stream, type, depth);
    stream << "},";

    for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
        print(stream, &*it, depth+1);
    }
}


void OutputJSONMachine::print_impl(ostream& stream, const AdHocAbstractRecord *type, unsigned int depth) {
	TypeBase::TypeHash hash=type->content_hash();
    if (doc_flags_.was_written(hash)) return;

    stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"AdHocAbstractRecord\"," << endl;
    stream << "\"parent\" : \"" << get_adhoc_parent_name(type) << "\"," << endl;
    type->write_attributes(stream);
    stream << "," << endl;

    print_abstract_record_keys(stream, dynamic_cast<const Type::AbstractRecord *>(type), depth);
    stream << "},";

    for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
        print(stream, &*it, depth+1);
    }

}



void OutputJSONMachine::print_abstract_record_keys(ostream& stream, const AbstractRecord *type, unsigned int depth) {

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





void OutputJSONMachine::print_impl(ostream& stream, const Selection *type, unsigned int depth) {
	TypeBase::TypeHash hash=type->content_hash();
    if (doc_flags_.was_written(hash)) return;

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
}


void OutputJSONMachine::print_impl(ostream& stream, const Integer *type, unsigned int depth) {
	TypeBase::TypeHash hash=type->content_hash();
    if (doc_flags_.was_written(hash)) return;

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
}


void OutputJSONMachine::print_impl(ostream& stream, const Double *type, unsigned int depth) {
	TypeBase::TypeHash hash=type->content_hash();
    if (doc_flags_.was_written(hash)) return;

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
}


void OutputJSONMachine::print_impl(ostream& stream, const Bool *type, unsigned int depth) {
	TypeBase::TypeHash hash=type->content_hash();
    if (doc_flags_.was_written(hash)) return;

    stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"Bool\"," << endl;
    stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	type->write_attributes(stream);
    stream << endl;
	stream << "},";
}


void OutputJSONMachine::print_impl(ostream& stream, const String *type, unsigned int depth) {
	TypeBase::TypeHash hash=type->content_hash();
    if (doc_flags_.was_written(hash)) return;

    stream << "{" << endl;
    stream << "\"id\" : \"" << format_hash(hash) << "\"," << endl;
    stream << "\"input_type\" : \"String\"," << endl;
    stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	type->write_attributes(stream);
    stream << endl;
	stream << "},";
}


void OutputJSONMachine::print_impl(ostream& stream, const FileName *type, unsigned int depth) {
	TypeBase::TypeHash hash=type->content_hash();
    if (doc_flags_.was_written(hash)) return;

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
}






std::ostream& operator<<(std::ostream& stream, OutputJSONMachine type_output) {
    return type_output.print(stream) << endl;
}


} // closing namespace Type
} // closing namespace Input
