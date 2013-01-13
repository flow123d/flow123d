/**
 * type_output.cc
 */

#include "input/type_output.hh"


namespace Input {
namespace Type {

using namespace std;


/*******************************************************************
 * implementation of OutputBase
 */

OutputBase::~OutputBase() {}


void OutputBase::print(ostream& stream) {
	doc_type_ = full_record;
	print(stream, type_);
}


void OutputBase::get_array_sizes(Array array, unsigned int &lower , unsigned int &upper ) {
	lower = array.data_->lower_bound_;
	upper = array.data_->upper_bound_;
}


void OutputBase::get_record_key(Record rec, unsigned int key_idx, Record::Key &key) {
	Record::KeyIter it = rec.begin() + key_idx;
	key = *it;
}


void OutputBase::get_integer_bounds(Integer integer, int &lower , int &upper ) {
	lower = integer.lower_bound_;
	upper = integer.upper_bound_;
}


void OutputBase::get_double_bounds(Double dbl, double &lower , double &upper ) {
	lower = dbl.lower_bound_;
	upper = dbl.upper_bound_;
}


void OutputBase::print(ostream& stream, const TypeBase *type, unsigned int depth) {
	if (typeid(*type) == typeid(Type::Record)) {
		print(stream, static_cast<const Type::Record *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::Array)) {
		print(stream, static_cast<const Type::Array *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::AbstractRecord)) {
		print(stream, static_cast<const Type::AbstractRecord *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::Selection)) {
		print(stream, static_cast<const Type::Selection *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::Integer)) {
		print(stream, static_cast<const Type::Integer *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::Double)) {
		print(stream, static_cast<const Type::Double *>(type), depth );
	} else
	if (typeid(*type) == typeid(Type::Bool)) {
		print(stream, static_cast<const Type::Bool *>(type), depth );
	} else {
		const Type::FileName * file_name_type = dynamic_cast<const Type::FileName *>(type);
        if (file_name_type != NULL ) {
        	print(stream, file_name_type, depth );
        	return;
        }

		const Type::String * string_type = dynamic_cast<const Type::String *>(type);
        if (string_type != NULL ) {
        	print(stream, string_type, depth );
        	return;
        }

        // default -> error
        xprintf(Err,"Unknown descendant of TypeBase class, name: %s\n", typeid(type).name());
	}
}


void OutputBase::write_value(std::ostream& stream, Default dft) {
	if (dft.is_obligatory() || dft.is_optional()) {
		stream << "<" << dft.value() << ">";
	} else {
		stream << "\"" << dft.value() << "\"";
	}
}




/*******************************************************************
 * implementation of OutputText
 */

void OutputText::print(ostream& stream, const Record *type, unsigned int depth) {
	if (! type->is_finished()) {
		xprintf(Warn, "Printing documentation of unfinished Input::Type::Record!\n");
	}

	switch (doc_type_) {
	case key_record:
		stream << "" << "Record '" << type->type_name() << "' (" << type->size() << " keys).";
		break;
	case full_record:
		if (! type->made_extensive_doc()) {
			type->set_made_extensive_doc(true);

			// header
			stream << endl;
			stream << "" << "Record '" << type->type_name() << "'";

			// reducible to key
			Record::KeyIter key_it = type->auto_conversion_key_iter();
			if (key_it != type->end()) {
				stream << ", reducible to key '" << key_it->key_ << "'";
			}

			stream << "" << " (" << type->size() << " keys).";
		    stream << endl;
		    stream << "" << "# " << type->description() << endl;
		    stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
		    // keys
		    doc_type_ = key_record;
		    for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
		    	size_setw_ = it->key_.size() + 3;
		        stream << setw(padding_size) << "" << it->key_ << " = ";
		        write_value(stream, it->default_);
		        stream << " is ";
		        print(stream, it->type_.get());
		        write_description(stream, it->description_);
		        stream << endl;

		    }
		    stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type->type_name() << endl;

		    // Full documentation of embedded record types.
		    doc_type_ = full_record;
		    if (depth_ == 0 || depth_ > depth) {
			    for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
			    	print(stream, it->type_.get(), depth+1);
			    }
		    }
		}
		break;
	}
}


void OutputText::print(ostream& stream, const Array *type, unsigned int depth) {

	switch (doc_type_) {
	case key_record:
		unsigned int lower_size, upper_size;

		get_array_sizes(*type, lower_size, upper_size);
		stream << "Array, size limits: [" << lower_size << ", " << upper_size << "] of type: " << endl;
		stream << setw(padding_size + size_setw_) << "";
		print(stream, type->data_->type_of_values_.get());
		break;
	case full_record:
		print(stream, type->data_->type_of_values_.get(), depth);
		break;
	}
}


void OutputText::print(ostream& stream, const AbstractRecord *type, unsigned int depth) {
	// Print documentation of abstract record
	switch (doc_type_) {
	case key_record:
		stream << "AbstractRecord '" << type->type_name() << "' with "<< type->child_size() << " descendants.";
		break;
	case full_record:
        if (! type->made_extensive_doc()) {

            // Extensive description
            type->set_made_extensive_doc(true);

            // header
            stream << endl;
            stream << "" << "AbstractRecord '" << type->type_name() << "' with " << type->child_size() << " descendants.";
            stream << endl;
            stream << "" << "# " << type->description() << endl;
            stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
            // descendants
            doc_type_ = key_record;
            for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
            	size_setw_ = 0;
                stream << setw(padding_size) << "";
                stream << "" << "Record '" << (*it).type_name() << "'";
                write_description(stream, it->description());
                stream << endl;
            }
            stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type->type_name() << endl;

            // Full documentation of embedded record types.
            doc_type_ = full_record;
            if (depth_ == 0 || depth_ > depth) {
                for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
                    print(stream, &*it, depth+1);
                }
            }
        }
		break;
	}
}


void OutputText::print(ostream& stream, const Selection *type, unsigned int depth) {
	if (! type->is_finished()) {
		xprintf(Warn, "Printing documentation of unfinished Input::Type::Selection!\n");
	}

	switch (doc_type_) {
	case key_record:
		stream << "Selection '" << type->type_name() << "' of " << type->size() << " values.";
		break;
	case full_record:
		if (! type->made_extensive_doc()) {
			type->set_made_extensive_doc(true);

			stream << endl << "Selection '" << type->type_name() << "' of " << type->size() << " values." << endl;
		    stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
		    // keys
		    for (Selection::keys_const_iterator it = type->begin(); it != type->end(); ++it) {
		        stream << setw(padding_size) << "" << it->key_ << " = " << it->value;
		        if (it->description_ != "") {
		        	stream << endl;
		        	stream << setw(2 * padding_size) << "" << "# " << it->description_ << "";
		        }
		        stream << endl;
		    }
		    stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type->type_name() << endl;
		}
		break;
	}
}


void OutputText::print(ostream& stream, const Integer *type, unsigned int depth) {
	if (doc_type_ == key_record) {
		int lower_bound, upper_bound;
		get_integer_bounds(*type, lower_bound, upper_bound);
		stream << "Integer in [" << lower_bound << ", " << upper_bound << "]";
	}
}


void OutputText::print(ostream& stream, const Double *type, unsigned int depth) {
	if (doc_type_ == key_record) {
		double lower_bound, upper_bound;
		get_double_bounds(*type, lower_bound, upper_bound);
		stream << "Double in [" << lower_bound << ", " << upper_bound << "]";
	}
}


void OutputText::print(ostream& stream, const Bool *type, unsigned int depth) {
	if (doc_type_ == key_record) {
		stream << "Bool";
	}
}


void OutputText::print(ostream& stream, const String *type, unsigned int depth) {
	if (doc_type_ == key_record) {
		stream << "String (generic)";
	}
}


void OutputText::print(ostream& stream, const FileName *type, unsigned int depth) {
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


void OutputText::write_description(std::ostream& stream, const string& str) {
    boost::tokenizer<boost::char_separator<char> > line_tokenizer(str, boost::char_separator<char>("\n"));
    boost::tokenizer<boost::char_separator<char> >::iterator it;

	// For every \n add padding at beginning of the next line.
	for(it = line_tokenizer.begin(); it != line_tokenizer.end(); ++it) {
		stream << endl;
		stream << setw(padding_size + size_setw_) << "" << "# " << *it;
	}
}


std::ostream& operator<<(std::ostream& stream, OutputText type_output) {
	type_output.print(stream);
	stream << endl;
	return stream;
}



/*******************************************************************
 * implementation of OutputJSONTemplate
 */


void OutputJSONTemplate::print(ostream& stream, const Record *type, unsigned int depth) {
	stream << endl;
	stream << setw(depth * padding_size) << "";
	stream << "# record " << type->type_name();
	if (type_name_.size()) {
		write_description(stream, description_);
		stream << endl << setw(depth * padding_size) << "" << type_name_ << " = ";
	} else {
		stream << endl << setw(depth * padding_size) << "";
	}

	stream << "{" << endl;
	for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
    	if (it->key_ == "TYPE") {
    		stream << setw((depth + 1) * padding_size) << "" << "TYPE = \"" << type->type_name() << "\"";
    	} else {
			type_name_ = it->key_;
			description_ = it->description_;
			size_setw_ = depth+1;
			value_ = it->default_;
			print(stream, it->type_.get(), depth+1);
		}
		stream << endl;
	}
	stream << setw(depth * padding_size) << "" << "}";
	if (depth == 0) {
		stream << endl;
	}
}


void OutputJSONTemplate::print(ostream& stream, const Array *type, unsigned int depth) {
	unsigned int lower_size, upper_size;
	get_array_sizes(*type, lower_size, upper_size);

	stream << endl;
	stream << setw(depth * padding_size) << "" << "# Array, size limits: [";
	stream << lower_size << ", " << upper_size << "] ";
	write_description(stream, description_);
	stream << endl;

	stream << setw(depth * padding_size) << "" << type_name_ << " = ";
	type_name_ = "";
	size_setw_ = depth + 1;
	stream << "[" << endl;

	print(stream, type->data_->type_of_values_.get(), depth+1);
	if (lower_size > 1) {
		stream << "," << endl;
		stream << setw((depth + 1) * padding_size) << "" << "< ";
		if (lower_size == 2) stream << "1 more entry >";
		else stream << (lower_size-1) << " more entries >";
	}

	stream << endl;
	stream << setw(depth * padding_size) << "" << "]";
}


void OutputJSONTemplate::print(ostream& stream, const AbstractRecord *type, unsigned int depth) {
	string rec_name = type_name_;

	stream << endl;
	stream << setw(depth * padding_size) << "" << "# abstract record " << type->type_name();
	write_description(stream, description_);
	stream << endl;
	stream << setw(depth * padding_size) << "";
	stream << "# " << std::setfill('-') << setw(20) << "" << std::setfill(' ') << " DESCENDANTS FOLLOWS";

    for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
    	type_name_ = rec_name;
    	description_ = it->description();
    	size_setw_ = depth;

    	if (it != type->begin_child_data()) {
	    	stream << ",";
		}
    	print(stream, &*it, depth);
    }

    stream << endl;
}


void OutputJSONTemplate::print(ostream& stream, const Selection *type, unsigned int depth) {
	stream << endl;
	stream << setw(depth * padding_size) << "" << "# Selection of " << type->size() << " values.";
	write_description(stream, description_);

	for (Selection::keys_const_iterator it = type->begin(); it != type->end(); ++it) {
		stream << endl;
		stream << setw(depth * padding_size) << "" << type_name_ << " = ";
		stream << "\"" << it->key_ << "\"";
        if (it->description_ != "") {
        	stream << setw(padding_size) << "" << "# " << it->description_ << "";
        }
    }

}


void OutputJSONTemplate::print(ostream& stream, const Integer *type, unsigned int depth) {
	// get bounds of integer
	int lower_bound, upper_bound;
	get_integer_bounds(*type, lower_bound, upper_bound);

	// test if value in value_.value() is not integer
	stringstream ss(value_.value());
	int i;
	bool invalid_val = (ss >> i).fail();

	// output
	stream << endl;
	stream << setw(depth * padding_size) << "" << "# Integer in [" << lower_bound << ", " << upper_bound << "]";
	write_description(stream, description_);
	stream << endl;
	stream << setw(depth * padding_size) << "";
	if (invalid_val) {
		stream << "# ";
	}
	if (type_name_.size()) {
		stream << "" << type_name_ << " = ";
	}
	if (invalid_val) {
		write_value(stream, value_);
	} else {
		stream << "" << value_.value() << "";
	}
}


void OutputJSONTemplate::print(ostream& stream, const Double *type, unsigned int depth) {
	// get bounds of double
	double lower_bound, upper_bound;
	get_double_bounds(*type, lower_bound, upper_bound);

	// test if value in value_.value() is not double
	stringstream ss(value_.value());
	double d;
	bool invalid_val = (ss >> d).fail();

	// output
	stream << endl;
	stream << setw(depth * padding_size) << "" << "# Double in [" << lower_bound << ", " << upper_bound << "]";
	write_description(stream, description_);
	stream << endl;
	stream << setw(depth * padding_size) << "";
	if (invalid_val) {
		stream << "# ";
	}
	if (type_name_.size()) {
		stream << "" << type_name_ << " = ";
	}
	if (invalid_val) {
		write_value(stream, value_);
	} else {
		stream << "" << value_.value() << "";
	}
}


void OutputJSONTemplate::print(ostream& stream, const Bool *type, unsigned int depth) {
	// test if in value_.value() is stored boolean value
	bool valid_val = (value_.value() == "true") || (value_.value() == "false");

	// output
	stream << endl;
	stream << setw(depth * padding_size) << "" << "# Boolean ";
	write_description(stream, description_);
	stream << endl;
	stream << setw(depth * padding_size) << "";
	if (!valid_val) {
		stream << "# ";
	}
	if (type_name_.size()) {
		stream << "" << type_name_ << " = ";
	}
	if (valid_val) {
		stream << "" << value_.value() << "";
	} else {
		write_value(stream, value_);
	}
}


void OutputJSONTemplate::print(ostream& stream, const String *type, unsigned int depth) {
	stream << endl;
	stream << setw(depth * padding_size) << "" << "# String ";
	write_description(stream, description_);
	stream << endl;
	stream << setw(depth * padding_size) << "";
	if (value_.is_obligatory() || value_.is_optional()) {
		stream << "# ";
	}
	if (type_name_.size()) {
		stream << "" << type_name_ << " = ";
	}
	write_value(stream, value_);
}


void OutputJSONTemplate::print(ostream& stream, const FileName *type, unsigned int depth) {
	stream << endl;
	stream << setw(depth * padding_size) << "" << "# FileName of ";

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

	write_description(stream, description_);
	stream << endl;
	stream << setw(depth * padding_size) << "";
	if (value_.is_obligatory() || value_.is_optional()) {
		stream << "# ";
	}
	if (type_name_.size()) {
		stream << "" << type_name_ << " = ";
	}
	write_value(stream, value_);
}


void OutputJSONTemplate::write_description(std::ostream& stream, const string& str) {
    boost::tokenizer<boost::char_separator<char> > line_tokenizer(str, boost::char_separator<char>("\n"));
    boost::tokenizer<boost::char_separator<char> >::iterator it;

	// For every \n add padding at beginning of the next line.
	for(it = line_tokenizer.begin(); it != line_tokenizer.end(); ++it) {
		stream << endl;
		stream << setw(size_setw_ * padding_size) << "" << "# " << *it;
	}
}


std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output) {
	type_output.print(stream);
	stream << endl;
	return stream;
}


} // closing namespace Type
} // closing namespace Input
