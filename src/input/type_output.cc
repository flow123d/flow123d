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
	doc_flags_.filter_ = NULL;
    TypeBase::lazy_finish();
}



ostream&  OutputBase::print(ostream& stream) {
	doc_type_ = full_record;
	doc_flags_.clear();

	print(stream, type_, depth_);
	return stream;
}


void OutputBase::set_filter(string regex_filter) {
    doc_flags_.reg_exp_ = regex_filter;
	doc_flags_.filter_ = new boost::regex(regex_filter);
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



void OutputBase::get_array_type(Array array, boost::shared_ptr<const TypeBase> &arr_type) {
    arr_type = array.data_->type_of_values_;
}



const string & OutputBase::get_record_description(const Record *rec) {
    return rec->data_->description_;
}



void OutputBase::get_record_key(Record rec, unsigned int key_idx, Record::Key &key) {
	Record::KeyIter it = rec.begin() + key_idx;
	key = *it;
}



void OutputBase::get_parent_ptr(Record rec, boost::shared_ptr<AbstractRecord> &parent_ptr) {
	parent_ptr = rec.data_->parent_ptr_;
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




const void * OutputBase::get_record_data(const Record *rec) {
	return rec->data_.get();
}

const void * OutputBase::get_abstract_record_data(const AbstractRecord *a_rec) {
	return a_rec->child_data_.get();
}

const void * OutputBase::get_selection_data(const Selection *sel) {
	return sel->data_.get();
}


const void * OutputBase::get_array_data(const Array *array) {
	return array->data_.get();
}


const void * OutputBase::get_type_base_data(const TypeBase *type) {
	if (typeid(*type) == typeid(Type::Record)) {
		return ( static_cast<const Type::Record *>(type) )->data_.get();
	} else
	if (typeid(*type) == typeid(Type::Array)) {
		return ( static_cast<const Type::Array *>(type) )->data_.get();
	} else
	if (typeid(*type) == typeid(Type::AbstractRecord)) {
		return ( static_cast<const Type::AbstractRecord *>(type) )->child_data_.get();
	} else
	if (typeid(*type) == typeid(Type::Selection)) {
		return ( static_cast<const Type::Selection *>(type) )->data_.get();
	}

	return NULL;
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
    boost::tokenizer<boost::char_separator<char> > line_tokenizer(str, boost::char_separator<char>("\n"));
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


OutputBase::ProcessedTypes::~ProcessedTypes() {
    if (filter_ != NULL) {
        delete filter_;
    }
}


void OutputBase::ProcessedTypes::clear() {
	key_to_index.erase(key_to_index.begin(), key_to_index.end());
	keys.erase(keys.begin(), keys.end());
	if (filter_ != NULL) {
		full_type_names.erase(full_type_names.begin(), full_type_names.end());
	}
}


unsigned int OutputBase::ProcessedTypes::type_index(const void * type_data) const {
    ProcessedTypes::key_to_index_const_iter it = key_to_index.find(type_data);
    if (it != key_to_index.end()) return it->second;

    //THROW( ExcRecordKeyNotFound() << EI_KeyName(key) << EI_Record(*this) );

    return keys.size();
}


/*void OutputBase::ProcessedTypes::set_reference(const void * type_data, const string& ref) {
	ProcessedTypes::key_to_index_const_iter data_it = key_to_index.find(type_data);

	ASSERT(data_it != key_to_index.end(), "Invalid key '%s' in OutputBase::OutputData object in set_reference method!\n", (static_cast<const Type::TypeBase *>(type_data))->type_name().c_str());

	KeyIter it = keys.begin() + data_it->second;
	(*it).reference_ = ref;
}*/


bool OutputBase::ProcessedTypes::was_written(const void * type_data, string full_name) {
	KeyIter it = keys.begin() + type_index(type_data);
	bool has_extensive = it != keys.end();

	if (filter_ == NULL) return has_extensive;

	std::string filtered = boost::regex_replace(full_name, *filter_, "");
	return (full_type_names.find(filtered) != full_type_names.end()) | has_extensive;
}


void OutputBase::ProcessedTypes::mark_written(const void *type_data, string full_name, string reference) {
	key_to_index_const_iter it = key_to_index.find(type_data);
	if ( it == key_to_index.end() ) {
	   key_to_index.insert( std::make_pair(type_data, keys.size()) );
	   Key tmp_key = { (unsigned int)keys.size(), type_data, reference};
	   keys.push_back(tmp_key);
	}

	if (filter_ != NULL) {
		std::string filtered = boost::regex_replace(full_name, *filter_, "");
		if ( full_type_names.find(filtered) == full_type_names.end() ) {
			full_type_names.insert(filtered);
		}
	}
}


const string OutputBase::ProcessedTypes::get_reference(const void * type_data) const {
	ProcessedTypes::key_to_index_const_iter data_it = key_to_index.find(type_data);

	ASSERT(data_it != key_to_index.end(), "Invalid key '%s' in OutputBase::OutputData object in get_reference method!\n", (static_cast<const Type::TypeBase *>(type_data))->type_name().c_str());

	KeyIter it = keys.begin()+data_it->second;
	return (*it).reference_;
}






/*******************************************************************
 * implementation of OutputText
 */

void OutputText::print_impl(ostream& stream, const Record *type, unsigned int depth) {
	if (! type->is_finished()) {
		xprintf(Warn, "Printing documentation of unfinished Input::Type::Record!\n");
	}

	switch (doc_type_) {
	case key_record:
		stream << "" << "Record '" << type->type_name() << "' (" << type->size() << " keys).";
		break;
	case full_record:
		const void * data_ptr = get_record_data(type);
		if (! doc_flags_.was_written(data_ptr, type->full_type_name())) {
			doc_flags_.mark_written(data_ptr, type->full_type_name());

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
		        print(stream, it->type_.get(), 0);
		        write_description(stream, it->description_, padding_size+size_setw_);
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


void OutputText::print_impl(ostream& stream, const Array *type, unsigned int depth) {
	boost::shared_ptr<const TypeBase> array_type;
	get_array_type(*type, array_type);

	switch (doc_type_) {
	case key_record:
		unsigned int lower_size, upper_size;

		get_array_sizes(*type, lower_size, upper_size);
		stream << "Array, size limits: [" << lower_size << ", " << upper_size << "] of type: " << endl;
		stream << setw(padding_size + size_setw_) << "" << "#### ";
		print(stream, array_type.get(), 0);
		break;
	case full_record:
		print(stream, array_type.get(), depth);
		break;
	}
}


void OutputText::print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth) {
	// Print documentation of abstract record
	switch (doc_type_) {
	case key_record:
		stream << "AbstractRecord '" << type->type_name() << "' with "<< type->child_size() << " descendants.";
		break;
	case full_record:
		const void * data_ptr = get_abstract_record_data(type);
		if (! doc_flags_.was_written(data_ptr, type->full_type_name()) ) {

            // Extensive description
			doc_flags_.mark_written(data_ptr, type->full_type_name());

            // header
            stream << endl;
            stream << "" << "AbstractRecord '" << type->type_name() << "' with " << type->child_size() << " descendants.";
            write_description(stream, OutputBase::get_record_description( type ), 0);
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
            if (depth_ == 0 || depth_ > depth) {
                for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
                    print(stream, &*it, depth+1);
                }
            }
        }
		break;
	}
}


void OutputText::print_impl(ostream& stream, const Selection *type, unsigned int depth) {
	if (! type->is_finished()) {
		xprintf(Warn, "Printing documentation of unfinished Input::Type::Selection!\n");
	}

	switch (doc_type_) {
	case key_record:
		stream << "Selection '" << type->type_name() << "' of " << type->size() << " values.";
		break;
	case full_record:
		const void * data_ptr = get_selection_data(type);

		if (! doc_flags_.was_written(data_ptr, type->full_type_name()) ) {
			doc_flags_.mark_written(data_ptr, type->full_type_name());

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


void OutputText::print_impl(ostream& stream, const Integer *type, unsigned int depth) {
	if (doc_type_ == key_record) {
		int lower_bound, upper_bound;
		get_integer_bounds(*type, lower_bound, upper_bound);
		stream << "Integer in [" << lower_bound << ", " << upper_bound << "]";
	}
}


void OutputText::print_impl(ostream& stream, const Double *type, unsigned int depth) {
	if (doc_type_ == key_record) {
		double lower_bound, upper_bound;
		get_double_bounds(*type, lower_bound, upper_bound);
		stream << "Double in [" << lower_bound << ", " << upper_bound << "]";
	}
}


void OutputText::print_impl(ostream& stream, const Bool *type, unsigned int depth) {
	if (doc_type_ == key_record) {
		stream << "Bool";
	}
}


void OutputText::print_impl(ostream& stream, const String *type, unsigned int depth) {
	if (doc_type_ == key_record) {
		stream << "String (generic)";
	}
}


void OutputText::print_impl(ostream& stream, const FileName *type, unsigned int depth) {
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
 * implementation of OutputJSONTemplate
 */

ostream& OutputJSONTemplate::print(ostream& stream) {
    key_name_ = "";
    return OutputBase::print(stream);
}


void OutputJSONTemplate::print_impl(ostream& stream, const Record *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record:
			stream << "# Record " << type->type_name();
			break;
		case full_record:
			stream << endl << setw(depth * padding_size) << "";
			if (key_name_.size()) {
				stream << key_name_ << " = ";
			}

			const void *data_ptr = get_record_data(type); // get pointer to type->data_
			if ( doc_flags_.was_written(data_ptr, type->full_type_name())) {
				stream << "{REF=\"" << doc_flags_.get_reference(data_ptr) << "\"}";
			} else {
				doc_flags_.mark_written( data_ptr, type->full_type_name(), reference_ );

				stream << "{";
				if (OutputBase::get_record_description(type).size()) {
					size_setw_ = depth+1;
					write_description(stream, OutputBase::get_record_description(type), padding_size*size_setw_, 2);
				}
				for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
					if (typeid(*(it->type_.get())) == typeid(Type::AbstractRecord)) {
						reference_ = doc_flags_.get_reference(data_ptr) + "/" + "#" + it->key_;
					} else if ( (typeid(*(it->type_.get())) == typeid(Type::Record))
							| (typeid(*(it->type_.get())) == typeid(Type::Array))
							| (typeid(*(it->type_.get())) == typeid(Type::Selection)) ) {
						reference_ = doc_flags_.get_reference(data_ptr) + "/" + it->key_;
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
						write_description(stream, it->description_, padding_size*size_setw_);
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


void OutputJSONTemplate::print_impl(ostream& stream, const Array *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record:
			unsigned int lower_size, upper_size;
			get_array_sizes(*type, lower_size, upper_size);

			stream << "# Array, size limits: [" << lower_size << ", " << upper_size << "] ";
			break;
		case full_record:
			bool has_opt_prefix = value_.is_optional() | value_.has_value_at_read_time(); // key contains OPT_ prefix
			boost::shared_ptr<const TypeBase> array_type;
			const void * data_ptr = get_array_data(type); // get pointer to type->data_

			get_array_type(*type, array_type);
			doc_flags_.mark_written( data_ptr, type->full_type_name(), reference_ );

			if ( (typeid(*(array_type.get())) == typeid(Type::Record))
					| (typeid(*(array_type.get())) == typeid(Type::AbstractRecord))
					| (typeid(*(array_type.get())) == typeid(Type::Array))) {
				reference_ = doc_flags_.get_reference(data_ptr) + "/0";
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
			print(stream, array_type.get(), depth+1);

			doc_type_ = full_record;
			if ( ! ( (typeid( *(array_type.get()) ) == typeid(Type::Integer))
					| (typeid( *(array_type.get()) ) == typeid(Type::Double))
					| (typeid( *(array_type.get()) ) == typeid(Type::Bool))
					| (typeid( *(array_type.get()) ) == typeid(Type::String))
					| (typeid( *(array_type.get()) ) == typeid(Type::FileName)) ) ) {
				print(stream, array_type.get(), depth+1);
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


void OutputJSONTemplate::print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record:
			stream << "# abstract record " << type->type_name();
			break;
		case full_record:
			string rec_name = key_name_;

			std::vector<string> refs;
			boost::split(refs, reference_, boost::is_any_of("#"));
		    ASSERT( refs.size() == 2, "Invalid reference of %s, size %d\n", type->type_name().c_str(), refs.size());

		    stream << endl;
			stream << setw(depth * padding_size) << "";
			stream << "# " << std::setfill('-') << setw(20) << "" << std::setfill(' ') << " DESCENDANTS FOLLOWS";

		    for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
		    	reference_ = refs[0] + it->type_name() + "_" + refs[1];

		    	key_name_ = it->type_name() + "_" + rec_name;
		    	size_setw_ = depth;

		    	if (it != type->begin_child_data()) {
			    	stream << ",";
				}
		    	doc_type_ = key_record;
		    	stream << endl;
		    	stream << setw((depth) * padding_size) << "";
		    	print(stream, &*it, depth);
		    	write_description(stream, OutputBase::get_record_description( &(*it) ), padding_size*size_setw_);
		    	doc_type_ = full_record;
		    	print(stream, &*it, depth);
		    }

			break;
	}

}


void OutputJSONTemplate::print_impl(ostream& stream, const Selection *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record: {
			unsigned int max_size = 0; // maximal size for setw of description

			stream << "# Selection of " << type->size() << " values";

			if (OutputBase::get_selection_description(type).size()) {
				//size_setw_ = depth+1;
				write_description(stream, OutputBase::get_selection_description(type), padding_size*depth, 2);
			}

			stream << endl << setw(depth * padding_size) << "" << "# Possible values:";

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


void OutputJSONTemplate::print_impl(ostream& stream, const Integer *type, unsigned int depth) {
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


void OutputJSONTemplate::print_impl(ostream& stream, const Double *type, unsigned int depth) {
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


void OutputJSONTemplate::print_impl(ostream& stream, const Bool *type, unsigned int depth) {
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


void OutputJSONTemplate::print_impl(ostream& stream, const String *type, unsigned int depth) {
	switch (doc_type_) {
		case key_record:
			stream << "# String ";
			break;
		case full_record:
			print_default_value(stream, depth, "\"\"", false, true);
			break;
	}

}


void OutputJSONTemplate::print_impl(ostream& stream, const FileName *type, unsigned int depth) {
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
		write_default_value(stream, value_);
	} else {
		stream << "" << value_.value() << "";
	}
}



/*******************************************************************
 * implementation of OutputLatex
 */

namespace internal {
class output_filter : public boost::iostreams::multichar_output_filter {
public:
    template<typename Sink>
    std::streamsize write(Sink& snk, const char* s, streamsize n)
    {
        std::streamsize n_out = 0;
        while (n != 0) {
            --n;
            if (s[0] == '_' || s[0] == '$') {
                boost::iostreams::put(snk,'\\');
            }
            boost::iostreams::put(snk, *s++); ++n_out;
        }
        return n_out;
    }
};


/**
 * Prints range specification for Array size, Integer, and Double. Omit natural bounds (very long),
 * Omit whole specification if  both limits are natural.
 */
template <class T>
ostream & print_range(ostream& stream, T range_min, T range_max) {
    T min_val = std::numeric_limits<T>::min();
    T max_val = std::numeric_limits<T>::max();
    if (min_val > 0) min_val = -max_val;

    if (range_min != min_val) {
        if (range_max != max_val) stream << "[" << range_min << ", " << range_max << "]";
        else stream << "[" << range_min << ", ]";
    } else {
        if (range_max != max_val) {
            cout << "DBG" << range_max << " " << max_val << " " << max_val-range_max << endl;
            stream << "[ ," << range_max << "]";
        }
    }
    return stream;
}

/**
 * Make \HTRaised{prefix::str}. Hyper target raised to scroll to correct position.
 */
std::string hyper_target( const std::string &prefix, const std::string &str) {
    string label=prefix + "::" + str;
    boost::replace_all(label, "_", "-");
    boost::replace_all(label, ">", "");
    // \hyperlink{<prefix>::str}{str}
    return "\\HTRaised{" + label + "}{" + str + "}";
}

/**
 * Make \hyperlink{prefix::str, str}.
 */
std::string hyper_link( const std::string &prefix, const std::string &str) {
    string label=prefix + "::" + str;
    boost::replace_all(label, "_", "-");
    boost::replace_all(label, ">", "");
    // \hyperlink{<prefix>::str}{str}
    return "\\hyperlink{" + label + "}{" + str +"}";
}

/**
 * Make bidirectional link.
 */
std::string hyper_B( const std::string &prefix, const std::string &str) {
    string label=prefix + "::" + str;
    boost::replace_all(label, "_", "-");
    boost::replace_all(label, ">", "");
    // \hyperlink{<prefix>::str}{str}
    return "\\hyperB{" + label + "}{" + str +"}";
}

} // namespace internal



ostream& OutputLatex::print(ostream& stream) {
    boost::iostreams::filtering_ostream out;
    out.push(internal::output_filter());
    out.push(stream);
    OutputBase::print(out);
    return stream;
}


void OutputLatex::print_impl(ostream& stream, const Record *type, unsigned int depth) {
    if (! type->is_finished()) {
        xprintf(Warn, "Printing documentation of unfinished Input::Type::Record!\n");
    }

    switch (doc_type_) {
    case key_record:
        stream << "record: " << internal::hyper_link("IT", type->type_name());
        break;
    case full_record:
    	const void * data_ptr = get_record_data(type);
    	if (! doc_flags_.was_written(data_ptr, type->full_type_name()) ) {
			doc_flags_.mark_written(data_ptr, type->full_type_name());

            // header
            stream << endl <<"\\begin{RecordType}{"
                   << internal::hyper_target("IT", type->type_name()) << "}";

            // parent record
            boost::shared_ptr<AbstractRecord> parent_ptr;
            get_parent_ptr(*type, parent_ptr);
            if (parent_ptr) {
                stream << "{" << internal::hyper_link("IT", parent_ptr->type_name()) <<"}";
            } else {
                stream << "{}";
            }

            // reducible to key
            Record::KeyIter key_it = type->auto_conversion_key_iter();
            if (key_it != type->end()) {
                stream << "{" << internal::hyper_link( type->type_name(), key_it->key_) << "}";
            } else {
                stream << "{}";
            }
            // add info and description
            stream << "{\\AddDoc{" << type->type_name() +"}}{"  << OutputBase::get_record_description(type) << "}";
            stream << endl;

            // keys
            doc_type_ = key_record;
            for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {

                stream << "\\KeyItem{" << internal::hyper_B( type->type_name(), it->key_) << "}";
                stream << "{";
                print(stream, it->type_.get(), 0);
                stream << "}";

                if (it->default_.is_obligatory()) {
                    stream << "{\\textless\\it obligatory\\textgreater}";
                } else if (it->default_.is_optional()) {
                    stream << "{\\textless\\it optional\\textgreater}";
                } else if (it->default_.has_value_at_read_time()) {
                    stream << "{\"" << it->default_.value() << "\"}";
                } else {
                    stream << "{" << it->default_.value() << "}";
                }
                
                string temp_desc = it->description_;
                boost::replace_all(temp_desc, "\n", "\\\\");
                stream << "{\\AddDoc{" << type->type_name() << "::" << it->key_ << "}}{"
                       << temp_desc << "}" << endl;
            }

            stream << "\\end{RecordType}" << endl;

            // Full documentation of embedded record types.
            doc_type_ = full_record;
            if (depth_ == 0 || depth_ > depth) {
                for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
                    if (it->key_ == "TYPE") continue;
                    print(stream, it->type_.get(), depth+1);
                }
            }
        }
        break;
    }
}




void OutputLatex::print_impl(ostream& stream, const Array *type, unsigned int depth) {
	boost::shared_ptr<const TypeBase> array_type;
	get_array_type(*type, array_type);

    switch (doc_type_) {
    case key_record:
        unsigned int lower_size, upper_size;

        get_array_sizes(*type, lower_size, upper_size);
        stream << "Array ";
        internal::print_range<unsigned int>(stream, lower_size, upper_size);
        stream << " of ";
        print(stream, array_type.get(), 0);
        break;
    case full_record:
        print(stream, array_type.get(), depth);
        break;
    }
}


void OutputLatex::print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth) {
    // Print documentation of abstract record
    switch (doc_type_) {
    case key_record:
        stream << "abstract type: " << internal::hyper_link("IT",type->type_name());
        break;
    case full_record:
    	const void * data_ptr = get_abstract_record_data(type);
    	if (! doc_flags_.was_written(data_ptr, type->full_type_name()) ) {

            // Extensive description
			doc_flags_.mark_written(data_ptr, type->full_type_name());

            // header
            stream << endl << "\\begin{AbstractType}{"
                   << internal::hyper_target("IT", type->type_name() ) << "}";
            const Record *default_desc = type->get_default_descendant();
            if (default_desc) {
                stream << "{" << internal::hyper_link( "IT", default_desc->type_name()) << "}";
            } else {
                stream << "{}";
            }
            // add info and description
            stream << "{\\AddDoc{" << type->type_name() << "}}{"  << OutputBase::get_record_description(type) << "}" << endl;

            // descendants
            doc_type_ = key_record;
            for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
                stream << "\\Descendant{" << internal::hyper_link( "IT", (*it).type_name() ) << "}" << endl;
            }
            stream << "\\end{AbstractType}" << endl;


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


void OutputLatex::print_impl(ostream& stream, const Selection *type, unsigned int depth) {
    if (! type->is_finished()) {
        xprintf(Warn, "Printing documentation of unfinished Input::Type::Selection!\n");
    }

    switch (doc_type_) {
    case key_record:
        if ( type->type_name().find("TYPE") != string::npos ) {
            stream<< "selection: " << type->type_name();
        } else {
            stream<< "selection: " << internal::hyper_link( "IT", type->type_name() );
        }
        break;
    case full_record:
    	const void * data_ptr = get_selection_data(type);
    	if (! doc_flags_.was_written(data_ptr, type->full_type_name()) ) {
			doc_flags_.mark_written(data_ptr, type->full_type_name());

            stream <<endl << "\\begin{SelectionType}{" << internal::hyper_target("IT", type->type_name() ) << "}";
            stream << "{" << OutputBase::get_selection_description(type) << "}" <<endl;
            // keys
            for (Selection::keys_const_iterator it = type->begin(); it != type->end(); ++it) {
                stream << "\\KeyItem{" <<  ( it->key_ ) << "}{" << it->description_ << "}" << endl;
            }
            stream << "\\end{SelectionType}" << endl;
        }
        break;
    }
}


void OutputLatex::print_impl(ostream& stream, const Integer *type, unsigned int depth) {
    if (doc_type_ == key_record) {
        int lower_bound, upper_bound;
        get_integer_bounds(*type, lower_bound, upper_bound);
        stream << "Integer ";
        internal::print_range<int>(stream, lower_bound, upper_bound);
    }
}


void OutputLatex::print_impl(ostream& stream, const Double *type, unsigned int depth) {
    if (doc_type_ == key_record) {
        double lower_bound, upper_bound;
        get_double_bounds(*type, lower_bound, upper_bound);
        stream << "Double ";
        internal::print_range<double>(stream, lower_bound, upper_bound);
    }
}


void OutputLatex::print_impl(ostream& stream, const Bool *type, unsigned int depth) {
    if (doc_type_ == key_record) {
        stream << "Bool";
    }
}


void OutputLatex::print_impl(ostream& stream, const String *type, unsigned int depth) {
    if (doc_type_ == key_record) {
        stream << "String (generic)";
    }
}


void OutputLatex::print_impl(ostream& stream, const FileName *type, unsigned int depth) {
    if (doc_type_ == key_record) {
        switch (type->get_file_type()) {
        case ::FilePath::input_file:
            stream << "input file name";
            break;
        case ::FilePath::output_file:
            stream << "output file name";
            break;
        }
    }
}



/*******************************************************************
 * implementation of OutputJSONMachine
 */

void OutputJSONMachine::print_impl(ostream& stream, const Record *type, unsigned int depth) {
	stream << "{" << endl;
	stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	stream << "\"full_name\" : \"" << type->full_type_name() << "\"," << endl;
	stream << "\"type\" : \"Record\"," << endl;
	stream << "\"description\" : \"" << boost::regex_replace(OutputBase::get_record_description(type), boost::regex("\\n"), "\\\\n") << "\"," << endl;

	// parent record
	boost::shared_ptr<AbstractRecord> parent_ptr;
	get_parent_ptr(*type, parent_ptr);
	stream << "\"parent\" : \"";
	if (parent_ptr) {
		stream << parent_ptr->type_name();
	}
	stream << "\"," << endl;

	// reducible to key
	Record::KeyIter key_it = type->auto_conversion_key_iter();
	stream << "\"reducible_to_key\" : \"";
	if (key_it != type->end()) {
		stream << key_it->key_;
	}
	stream << "\"," << endl;

	stream << "\"keys\" : [" << endl;

	for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
		string dft_type, dft_value;
		get_default(it, dft_type, dft_value);

		if (it != type->begin()) {
			stream << "," << endl;
		}
		stream << "{ \"key\" : \"" << it->key_ << "\"," << endl;
		stream << "\"description\" : \"" << boost::regex_replace(it->description_, boost::regex("\\n"), "\\\\n") << "\"," << endl;
		stream << "\"default\" : { \"type\" : \"" << dft_type << "\",  \"value\" : \"" << dft_value << "\" }," << endl;
		stream << "\"type\" : ";
		print(stream, it->type_.get(), 0);
		stream << "}";
	}

	stream << "]" << endl;
	stream << "}";
}


void OutputJSONMachine::print_impl(ostream& stream, const Array *type, unsigned int depth) {
    unsigned int lower_size, upper_size;
	boost::shared_ptr<const TypeBase> array_type;

    get_array_sizes(*type, lower_size, upper_size);
	get_array_type(*type, array_type);

	stream << "{" << endl;
	stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	stream << "\"full_name\" : \"" << type->full_type_name() << "\"," << endl;
	stream << "\"type\" : \"Array\"," << endl;
	stream << "\"range\" : [" << lower_size << ", " << upper_size << "]," << endl;
	stream << "\"subtype\" : ";
	print(stream, array_type.get(), 0);
	stream << endl << "}";
}


void OutputJSONMachine::print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth) {
	stream << "{" << endl;
	stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	stream << "\"full_name\" : \"" << type->full_type_name() << "\"," << endl;
	stream << "\"type\" : \"AbstractRecord\"," << endl;
	stream << "\"description\" : \"" << boost::regex_replace( OutputBase::get_record_description(type), boost::regex("\\n"), "\\\\n") << "\"," << endl;

	// default descendant
	const Record * desc = type->get_default_descendant();
	stream << "\"default_descendant\" : \"";
	if (desc) {
		stream << desc->type_name();
	}
	stream << "\"," << endl;

	stream << "\"keys\" : [" << endl;

	for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
		if (it != type->begin_child_data()) {
			stream << "," << endl;
		}

		stream << "{ \"key\" : \"" << it->type_name() << "\"," << endl;
		stream << "\"description\" : \"" << boost::regex_replace( OutputBase::get_record_description( &(*it) ), boost::regex("\\n"), "\\\\n") << "\"," << endl;
		stream << "\"type\" : ";
		print(stream, &*it, depth);
		stream << "}";
	}

	stream << "]" << endl;
	stream << "}";
}


void OutputJSONMachine::print_impl(ostream& stream, const Selection *type, unsigned int depth) {
	stream << "{" << endl;
	stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	stream << "\"full_name\" : \"" << type->full_type_name() << "\"," << endl;
	stream << "\"type\" : \"Selection\"," << endl;
	stream << "\"description\" : \"" << boost::regex_replace( OutputBase::get_selection_description(type), boost::regex("\\n"), "\\\\n") << "\"," << endl;

	stream << "\"values\" : [" << endl;

	for (Selection::keys_const_iterator it = type->begin(); it != type->end(); ++it) {
		if (it != type->begin()) {
			stream << "," << endl;
		}
		stream << "{ \"value\" : \"" << it->value << "\", \"name\" : \"" << it->key_ << "\", \"description\" : \"" << it->description_ << "\" }";
	}

	stream << "]" << endl;
	stream << "}";
}


void OutputJSONMachine::print_impl(ostream& stream, const Integer *type, unsigned int depth) {
	int lower, upper;
	get_integer_bounds(*type, lower, upper);

	stream << "{" << endl;
	stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	stream << "\"full_name\" : \"" << type->full_type_name() << "\"," << endl;
	stream << "\"type\" : \"Integer\"," << endl;
	stream << "\"range\" : [" << lower << ", " << upper << "]" << endl;
	stream << "}";
}


void OutputJSONMachine::print_impl(ostream& stream, const Double *type, unsigned int depth) {
	double lower, upper;
	get_double_bounds(*type, lower, upper);

	stream << "{" << endl;
	stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	stream << "\"full_name\" : \"" << type->full_type_name() << "\"," << endl;
	stream << "\"type\" : \"Double\"," << endl;
	stream << "\"range\" : [" << lower << ", " << upper << "]" << endl;
	stream << "}";
}


void OutputJSONMachine::print_impl(ostream& stream, const Bool *type, unsigned int depth) {
	stream << "{" << endl;
	stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	stream << "\"full_name\" : \"" << type->full_type_name() << "\"," << endl;
	stream << "\"type\" : \"Bool\"" << endl;
	stream << "}";
}


void OutputJSONMachine::print_impl(ostream& stream, const String *type, unsigned int depth) {
	stream << "{" << endl;
	stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	stream << "\"full_name\" : \"" << type->full_type_name() << "\"," << endl;
	stream << "\"type\" : \"String\"" << endl;
	stream << "}";
}


void OutputJSONMachine::print_impl(ostream& stream, const FileName *type, unsigned int depth) {
	stream << "{" << endl;
	stream << "\"name\" : \"" << type->type_name() << "\"," << endl;
	stream << "\"full_name\" : \"" << type->full_type_name() << "\"," << endl;

	stream << "\"type\" : \"FileName of ";

	switch (type->get_file_type()) {
	case ::FilePath::input_file:
		stream << "input file\"";
		break;
	case ::FilePath::output_file:
		stream << "output file\"";
		break;
	default:
		stream << "file with unknown type\"";
		break;
	}

	stream << "}";
}




std::ostream& operator<<(std::ostream& stream, OutputText type_output) {
    return type_output.print(stream) << endl;
}



std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output) {
    return type_output.print(stream) << endl;
}



std::ostream& operator<<(std::ostream& stream, OutputLatex type_output) {
    return type_output.print(stream) << endl;
}

std::ostream& operator<<(std::ostream& stream, OutputJSONMachine type_output) {
    return type_output.print(stream) << endl;
}


} // closing namespace Type
} // closing namespace Input
