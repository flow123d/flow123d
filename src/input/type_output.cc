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


void OutputBase::write_description(std::ostream& stream, const string& str) {
    boost::tokenizer<boost::char_separator<char> > line_tokenizer(str, boost::char_separator<char>("\n"));
    boost::tokenizer<boost::char_separator<char> >::iterator it;

	// For every \n add padding at beginning of the next line.
	for(it = line_tokenizer.begin(); it != line_tokenizer.end(); ++it) {
		stream << setw(padding_size) << "" << "# " << *it << endl;
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
		stream << "" << "Record '" << type->type_name() << "' with " << type->size() << " keys.";
		break;
	case full_record:
		// header
		stream << endl;
	    stream << "" << "Record '" << type->type_name() << "' with " << type->size() << " keys.";
	    stream << endl;
	    stream << "" << "# " << type->description() << endl;
	    stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
	    // keys
	    doc_type_ = key_record;
	    for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
	        stream << setw(padding_size) << "" << it->key_ << " = <" << it->default_.value() << "> is ";
	        print(stream, it->type_.get());
	        stream << endl;
	        write_description(stream, it->description_);

	    }
	    stream << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type->type_name() << endl;

	    // Full documentation of embedded record types.
	    doc_type_ = full_record;
	    if (depth_ == 0 || depth_ > depth) {
		    for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
		    	print(stream, it->type_.get(), depth+1);
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
		stream << setw(padding_size) << "";
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
                stream << setw(padding_size) << "";
            	print(stream, &*it);
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
		        if (it->description_ != "")
		            stream << " (" << it->description_ << ")";
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


std::ostream& operator<<(std::ostream& stream, OutputText type_output) {
	type_output.print(stream);
	stream << endl;
	return stream;
}



/*******************************************************************
 * implementation of OutputJSONTemplate
 */


void OutputJSONTemplate::print(ostream& stream, const Record *type, unsigned int depth) {
	stream << "{" << endl;
	for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
		stream << setw((depth+1) * padding_size) << "" << it->key_ << " = ";
		print(stream, it->type_.get(), depth+1);
		stream << endl;
	}
	stream << setw(depth * padding_size) << "" << "}";
	if (depth == 0) {
		stream << endl;
	}
}


void OutputJSONTemplate::print(ostream& stream, const Array *type, unsigned int depth) {
	unsigned int lower_size, upper_size, minimum=2;
	get_array_sizes(*type, lower_size, upper_size);
	lower_size = std::max(lower_size, minimum);

	stream << "[" << endl;

	for (unsigned int i=0; i<lower_size; i++) {
		if (i > 0) {
			stream << "," << endl;
		}
		stream << setw((depth + 1) * padding_size) << "";
		print(stream, type->data_->type_of_values_.get(), depth+1);
	}

	stream << endl;
	stream << setw(depth * padding_size) << "" << "]";
}


void OutputJSONTemplate::print(ostream& stream, const AbstractRecord *type, unsigned int depth) {
	//aaa
}


void OutputJSONTemplate::print(ostream& stream, const Selection *type, unsigned int depth) {
	stream << "<";
	for (Selection::keys_const_iterator it = type->begin(); it != type->end(); ++it) {
		if (it != type->begin()) {
			stream << " | ";
		}
		stream << "\"" << it->key_ << "\"";
    }
	stream << ">";

}


void OutputJSONTemplate::print(ostream& stream, const Integer *type, unsigned int depth) {
	stream << "1";
}


void OutputJSONTemplate::print(ostream& stream, const Double *type, unsigned int depth) {
	stream << "2.5";
}


void OutputJSONTemplate::print(ostream& stream, const Bool *type, unsigned int depth) {
	stream << "false";
}


void OutputJSONTemplate::print(ostream& stream, const String *type, unsigned int depth) {
	stream << "\"Some string.\"";
}


void OutputJSONTemplate::print(ostream& stream, const FileName *type, unsigned int depth) {
	stream << "\"./mesh.msh\"";
}


std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output) {
	type_output.print(stream);
	stream << endl;
	return stream;
}


} // closing namespace Type
} // closing namespace Input
