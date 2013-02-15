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



OutputBase::OutputBase(const TypeBase *type, unsigned int depth)
: type_(type), depth_(depth)
{
    TypeBase::lazy_finish();
}



void OutputBase::print(ostream& stream) {
	doc_type_ = full_record;
	type_->reset_doc_flags();
	print(stream, type_, depth_);
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

			// parent record
			if (type->data_->parent_ptr_) {
				stream << ", implementation of " << type->data_->parent_ptr_->type_name();
			}

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
		        stream << endl;
		        stream << setw(padding_size + size_setw_) << "" <<"#### is ";
		        print(stream, it->type_.get(), 0);
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
		stream << setw(padding_size + size_setw_) << "" << "#### ";
		print(stream, type->data_->type_of_values_.get(), 0);
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
		        	stream << setw(padding_size + it->key_.size() + 3) << "" << "# " << it->description_ << "";
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


void OutputText::write_description(std::ostream& stream, const string& str, unsigned int hash_count) {
    boost::tokenizer<boost::char_separator<char> > line_tokenizer(str, boost::char_separator<char>("\n"));
    boost::tokenizer<boost::char_separator<char> >::iterator it;

	// For every \n add padding at beginning of the next line.
	for(it = line_tokenizer.begin(); it != line_tokenizer.end(); ++it) {
		stream << endl;
		stream << setw(padding_size + size_setw_) << "";
		stream << std::setfill('#') << setw(hash_count) << "" << std::setfill(' ') << " " << *it;
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
	switch (doc_type_) {
		case key_record:
			stream << "# Record " << type->type_name();
			break;
		case full_record:
			stream << endl << setw(depth * padding_size) << "";
			if (key_name_.size()) {
				stream << key_name_ << " = ";
			}
			if (type->made_extensive_doc()) {
				stream << "{REF=\"" << type->reference() << "\"}";
			} else {
				type->set_made_extensive_doc(true);
				stream << "{";
				if (type->description().size()) {
					size_setw_ = depth+1;
					write_description(stream, type->description(), 2);
				}
				for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
					if (typeid(*(it->type_.get())) == typeid(Type::AbstractRecord)) {
						it->type_.get()->set_reference( type->reference() + "#" + it->key_ + "/" );
					} else if ( !(typeid(*(it->type_.get())) == typeid(Type::Record))
							| (it->type_.get()->reference().size() == 0) ) {
						it->type_.get()->set_reference( type->reference() + it->key_ + "/" );
					}
					if (it != type->begin()) {
						stream << ",";
					}
					stream << endl;
			    	if (it->key_ == "TYPE") {
			    		stream << endl;
			    		stream << setw((depth + 1) * padding_size) << "" << "TYPE = \"" << type->type_name() << "\"";
			    	} else {
			    		key_name_ = it->key_;
						size_setw_ = depth+1;
						value_ = it->default_;

						doc_type_ = key_record;
						stream << endl;
						stream << setw((depth + 1) * padding_size) << "";
						print(stream, it->type_.get(), depth+1);
						write_description(stream, it->description_);
						doc_type_ = full_record;
						print(stream, it->type_.get(), depth+1);
					}
				}
				stream << endl;
				stream << setw(depth * padding_size) << "" << "}";
			}

			if (depth == 0) {
				stream << endl;
			}
			break;
	}

}


void OutputJSONTemplate::print(ostream& stream, const Array *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record:
			unsigned int lower_size, upper_size;
			get_array_sizes(*type, lower_size, upper_size);

			stream << "# Array, size limits: [" << lower_size << ", " << upper_size << "] ";
			break;
		case full_record:
			bool has_opt_prefix = value_.is_optional() | value_.has_value_at_read_time(); // key contains OPT_ prefix

			if (type->data_->type_of_values_.get()->reference().size() == 0) {
				type->data_->type_of_values_.get()->set_reference( type->reference() );
			}

			if (has_opt_prefix) {
				stream << endl;
				stream << setw(depth * padding_size) << "" << "# Optional key";
			}

			stream << endl;
			stream << setw(depth * padding_size) << "";

			if (key_name_.size()) {
				if (has_opt_prefix) {
					stream << "OPT_";
				}
				stream << key_name_ << " = " << "[";
			} else {
				stream << "[";
			}

			key_name_ = "";
			size_setw_ = depth + 1;

			doc_type_ = key_record;
			stream << endl;
			stream << setw((depth + 1) * padding_size) << "";
			print(stream, type->data_->type_of_values_.get(), depth+1);

			doc_type_ = full_record;
			if ( ! (typeid( *(type->data_->type_of_values_.get()) ) == typeid(Type::Integer)
					| typeid( *(type->data_->type_of_values_.get()) ) == typeid(Type::Double)
					| typeid( *(type->data_->type_of_values_.get()) ) == typeid(Type::Bool)
					| typeid( *(type->data_->type_of_values_.get()) ) == typeid(Type::String)
					| typeid( *(type->data_->type_of_values_.get()) ) == typeid(Type::FileName)) ) {
				print(stream, type->data_->type_of_values_.get(), depth+1);
			}
			//if (lower_size > 1) {
			//	stream << "," << endl;
			//	stream << setw((depth + 1) * padding_size) << "" << "< ";
			//	if (lower_size == 2) stream << "1 more entry >";
			//	else stream << (lower_size-1) << " more entries >";
			//}

			stream << endl;
			stream << setw(depth * padding_size) << "" << "]";
			break;
	}

}


void OutputJSONTemplate::print(ostream& stream, const AbstractRecord *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record:
			stream << "# abstract record " << type->type_name();
			break;
		case full_record:
			string rec_name = key_name_;

			std::vector<string> refs;
			string reference(type->reference());
			boost::split(refs, reference, boost::is_any_of("#"));
		    ASSERT( refs.size() == 2, "Invalid reference of %s, size %d\n", type->type_name().c_str(), refs.size());

			stream << endl;
			stream << setw(depth * padding_size) << "";
			stream << "# " << std::setfill('-') << setw(20) << "" << std::setfill(' ') << " DESCENDANTS FOLLOWS";

		    for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
		    	if (it->reference().size() == 0) {
		    		it->set_reference( refs[0] + it->type_name() + "_" + refs[1] );
		    	}

		    	key_name_ = it->type_name() + "_" + rec_name;
		    	size_setw_ = depth;

		    	if (it != type->begin_child_data()) {
			    	stream << ",";
				}
		    	doc_type_ = key_record;
		    	stream << endl;
		    	stream << setw((depth) * padding_size) << "";
		    	print(stream, &*it, depth);
		    	write_description(stream, it->description());
		    	doc_type_ = full_record;
		    	print(stream, &*it, depth);
		    }
		    type->set_reference( refs[0] + refs[1] );

		    stream << endl;
			break;
	}

}


void OutputJSONTemplate::print(ostream& stream, const Selection *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record: {
			unsigned int max_size = 0; // maximal size for setw of description

			stream << "# Selection of " << type->size() << " values:";

			for (Selection::keys_const_iterator it = type->begin(); it != type->end(); ++it) {
				max_size = std::max(max_size, (unsigned int)(it->key_.size()) );
			}

			for (Selection::keys_const_iterator it = type->begin(); it != type->end(); ++it) {
				stream << endl;
				stream << setw(depth * padding_size) << "" << "# \"" << it->key_ << "\"";
		        if (it->description_ != "") {
		        	stream << setw(max_size - it->key_.size()) << "" << " - " << it->description_ << "";
		        }
		    }

			stream << endl;
			stream << setw(depth * padding_size) << "";
			stream << "# " << std::setfill('-') << setw(10) << "" << std::setfill(' ');
			break;
		}
		case full_record:
			print_default_value(stream, depth, "\"\"", false, true);
			//if (value_.is_optional()) {
			//	stream << setw(depth * padding_size) << "" << "OPT_" << key_name_ << " = \"\"" ;
			//} else {
			//	stream << setw(depth * padding_size) << "" << key_name_ << " = \"" << value_.value()<< "\"" ;
			//}
			break;
	}

}


void OutputJSONTemplate::print(ostream& stream, const Integer *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record:
			int lower_bound, upper_bound;
			get_integer_bounds(*type, lower_bound, upper_bound);
			stream << "# Integer in [" << lower_bound << ", " << upper_bound << "]";
			break;
		case full_record:
			// test if value in value_.value() is not integer
			stringstream ss(value_.value());
			int i;
			bool invalid_val = (ss >> i).fail();

			print_default_value(stream, depth, "0", invalid_val);
			break;
	}

}


void OutputJSONTemplate::print(ostream& stream, const Double *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record:
			double lower_bound, upper_bound;
			get_double_bounds(*type, lower_bound, upper_bound);

			stream << "# Double in [" << lower_bound << ", " << upper_bound << "]";
			break;
		case full_record:
			// test if value in value_.value() is not double
			stringstream ss(value_.value());
			double d;
			bool invalid_val = (ss >> d).fail();

			print_default_value(stream, depth, "0", invalid_val);
			break;
	}

}


void OutputJSONTemplate::print(ostream& stream, const Bool *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record:
			stream << "# Boolean ";
			break;
		case full_record:
			// test if in value_.value() is stored boolean value
			bool invalid_val = (value_.value() != "true") & (value_.value() != "false");

			print_default_value(stream, depth, "false", invalid_val);
			break;
	}

}


void OutputJSONTemplate::print(ostream& stream, const String *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record:
			stream << "# String ";
			break;
		case full_record:
			print_default_value(stream, depth, "\"\"", false, true);
			break;
	}

}


void OutputJSONTemplate::print(ostream& stream, const FileName *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record:
			stream << "# FileName of ";

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
			break;
		case full_record:
			print_default_value(stream, depth, "\"\"", false, true);
			break;
	}

}


void OutputJSONTemplate::write_description(std::ostream& stream, const string& str, unsigned int hash_count) {
    boost::tokenizer<boost::char_separator<char> > line_tokenizer(str, boost::char_separator<char>("\n"));
    boost::tokenizer<boost::char_separator<char> >::iterator it;

	// For every \n add padding at beginning of the next line.
	for(it = line_tokenizer.begin(); it != line_tokenizer.end(); ++it) {
		stream << endl;
		stream << setw(size_setw_ * padding_size) << "";
		stream << std::setfill('#') << setw(hash_count) << "" << std::setfill(' ') << " " << *it;
	}
}


void OutputJSONTemplate::print_default_value(ostream& stream, unsigned int depth, string empty_val, bool invalid_val, bool has_quote) {
	stream << endl;
	stream << setw(depth * padding_size) << "";
	if (value_.is_optional() | value_.has_value_at_read_time()) {
		// optional and read time values have comment and key prefix OPT_
		if (value_.is_optional()) {
			stream << "# Optional key";
		} else {
			stream << "# Read time value - " << value_.value();
		}
		stream << endl;
		stream << setw(depth * padding_size) << "" << "OPT_";
	} else if (invalid_val & !value_.is_obligatory()) {
		// comment of non obligatory invalid values
		stream << "# ";
	}
	if (key_name_.size()) {
		stream << "" << key_name_ << " = ";
	}

	// printout of value
	if (value_.is_optional() | value_.has_value_at_read_time()) {
		stream << empty_val;
	} else if (invalid_val | has_quote) {
		write_value(stream, value_);
	} else {
		stream << "" << value_.value() << "";
	}
}


std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output) {
	type_output.print(stream);
	stream << endl;
	return stream;
}


} // closing namespace Type
} // closing namespace Input
