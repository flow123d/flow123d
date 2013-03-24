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
	type_->reset_doc_flags();
	doc_flags_.clear();
	doc_flags_.declare_key(type_->type_name(), false, "");
	print(stream, type_, depth_);
	return stream;
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


void OutputBase::get_parent_ptr(Record rec, boost::shared_ptr<AbstractRecord> &parent_ptr) {
	parent_ptr = rec.data_->parent_ptr_;
}


void OutputBase::get_array_type(Array array, boost::shared_ptr<const TypeBase> &arr_type) {
	arr_type = array.data_->type_of_values_;
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


void OutputBase::write_value(std::ostream& stream, Default dft) {
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


inline unsigned int OutputBase::key_index(const string& key) const {
    KeyHash key_h = key_hash(key);
    OutputData::key_to_index_const_iter it = doc_flags_.key_to_index.find(key_h);
    if (it != doc_flags_.key_to_index.end()) return it->second;

    //THROW( ExcRecordKeyNotFound() << EI_KeyName(key) << EI_Record(*this) );

    return size();
}


inline OutputBase::KeyIter OutputBase::key_iterator(const string& key) const {
	return begin() + key_index(key);
}


inline bool OutputBase::has_key_iterator(const string& key, KeyIter &it) const {
    KeyHash key_h = key_hash(key);
    OutputData::key_to_index_const_iter data_it = doc_flags_.key_to_index.find(key_h);
    if (data_it == doc_flags_.key_to_index.end()) {
        return false;
    } else {
        it = begin()+data_it->second;
        return true;
    }
}


inline OutputBase::KeyIter OutputBase::begin() const {
	return doc_flags_.keys.begin();
}


inline OutputBase::KeyIter OutputBase::end() const {
	return doc_flags_.keys.end();
}


inline bool OutputBase::has_key(const string& key) const {
	return key_iterator(key) != end();
}


inline unsigned int OutputBase::size() const {
	return doc_flags_.keys.size();
}


inline bool OutputBase::has_key_extensive(const string& key) const {
	if ( has_key(key) ) {
		return ( *key_iterator(key) ).extensive_doc_;
	}
	return false;
}


inline const string OutputBase::get_reference(const string& key) const {
	KeyHash key_h = key_hash(key);
	OutputData::key_to_index_const_iter data_it = doc_flags_.key_to_index.find(key_h);

	ASSERT(data_it != doc_flags_.key_to_index.end(), "Invalid key '%s' in OutputBase::OutputData object in get_reference method!\n", key.c_str());

	KeyIter it = begin()+data_it->second;
	return (*it).reference_;
}


inline void OutputBase::set_reference(const string& key, const string& ref) {
	KeyHash key_h = key_hash(key);
	OutputData::key_to_index_const_iter data_it = doc_flags_.key_to_index.find(key_h);

	ASSERT(data_it != doc_flags_.key_to_index.end(), "Invalid key '%s' in OutputBase::OutputData object in set_reference method!\n", key.c_str());

	KeyIter it = begin()+data_it->second;
	(*it).reference_ = ref;
}


inline void OutputBase::set_extensive_flag(const string& key, bool val) {
	ASSERT(has_key(key), "Invalid key '%s' in OutputBase::OutputData object in set_extensive_flag method!\n", key.c_str());

	( *key_iterator(key) ).extensive_doc_ = true;
}


inline void OutputBase::remove_key(const string& key) {
	ASSERT(has_key(key), "Invalid key '%s' in OutputBase::OutputData object during removal the key!\n", key.c_str());

	doc_flags_.key_to_index.erase(doc_flags_.key_to_index.find(key));
	doc_flags_.keys.erase(doc_flags_.keys.begin() + key_index(key));
}




/*******************************************************************
 * implementation of OutputBase::OutputData
 */


bool OutputBase::OutputData::declare_key(string key, bool extensive_doc, string reference) {
	KeyHash key_h = key_hash(key);
	key_to_index_const_iter it = key_to_index.find(key_h);
	if ( it == key_to_index.end() ) {
	   key_to_index.insert( std::make_pair(key_h, keys.size()) );
	   Key tmp_key = { (unsigned int)keys.size(), key, extensive_doc, reference};
	   keys.push_back(tmp_key);
	   return true;
	}
	return false;
}


bool OutputBase::OutputData::declare_key(string key) {
	return declare_key(key, true, "");
}


void OutputBase::OutputData::clear() {
	key_to_index.erase(key_to_index.begin(), key_to_index.end());
	keys.erase(keys.begin(), keys.end());
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
		if (! type->made_extensive_doc()) {
			type->set_made_extensive_doc(true);

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
                write_description(stream, it->description(), padding_size+size_setw_);
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
    //type_->set_reference( "" );
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
			//if (type->made_extensive_doc()) {
			if ( has_key_extensive(type->type_name()) ) {
				//stream << "{REF=\"" << type->reference() << "\"}";
				stream << "{REF=\"" << get_reference(type->type_name()) << "\"}";
			} else {
				type->set_made_extensive_doc(true); //delete
				set_extensive_flag(type->type_name());
				string ref = get_reference(type->type_name());

				stream << "{";
				if (type->description().size()) {
					size_setw_ = depth+1;
					write_description(stream, type->description(), padding_size*size_setw_, 2);
				}
				for (Record::KeyIter it = type->begin(); it != type->end(); ++it) {
					/*if (typeid(*(it->type_.get())) == typeid(Type::AbstractRecord)) {
						it->type_.get()->set_reference( type->reference() + "/" + "#" + it->key_ );
					} else if ( !(typeid(*(it->type_.get())) == typeid(Type::Record))
							| (it->type_.get()->reference().size() == 0) ) {
						it->type_.get()->set_reference( type->reference() + "/" + it->key_ );
					} */
					if (typeid(*(it->type_.get())) == typeid(Type::AbstractRecord)) {
						doc_flags_.declare_key( it->type_.get()->type_name(), false, ref + "/" + "#" + it->key_ );
					} else if ( (typeid(*(it->type_.get())) == typeid(Type::Record))
							| (typeid(*(it->type_.get())) == typeid(Type::Array))
							| (typeid(*(it->type_.get())) == typeid(Type::Selection)) ) {
						doc_flags_.declare_key( it->type_.get()->type_name(), false, ref + "/" + it->key_ );
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

			get_array_type(*type, array_type);

			/*if (array_type.get()->reference().size() == 0) {
				array_type.get()->set_reference( type->reference() + "/0" );
			} */
			if ( (typeid(*(array_type.get())) == typeid(Type::Record))
					| (typeid(*(array_type.get())) == typeid(Type::AbstractRecord))
					| (typeid(*(array_type.get())) == typeid(Type::Array))) {
				doc_flags_.declare_key( array_type.get()->type_name(), false, get_reference(type->type_name()) + "/0");
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
			if ( ! (typeid( *(array_type.get()) ) == typeid(Type::Integer)
					| typeid( *(array_type.get()) ) == typeid(Type::Double)
					| typeid( *(array_type.get()) ) == typeid(Type::Bool)
					| typeid( *(array_type.get()) ) == typeid(Type::String)
					| typeid( *(array_type.get()) ) == typeid(Type::FileName)) ) {
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
			//string reference(type->reference());
			string reference(get_reference(type->type_name()));
			boost::split(refs, reference, boost::is_any_of("#"));
		    ASSERT( refs.size() == 2, "Invalid reference of %s, size %d\n", type->type_name().c_str(), refs.size());

		    stream << endl;
			stream << setw(depth * padding_size) << "";
			stream << "# " << std::setfill('-') << setw(20) << "" << std::setfill(' ') << " DESCENDANTS FOLLOWS";

		    for (AbstractRecord::ChildDataIter it = type->begin_child_data(); it != type->end_child_data(); ++it) {
		    	/*if (it->reference().size() == 0) {
		    		it->set_reference( refs[0] + it->type_name() + "_" + refs[1] );
		    	} */
		    	doc_flags_.declare_key( it->type_name(), false, refs[0] + it->type_name() + "_" + refs[1] );

		    	key_name_ = it->type_name() + "_" + rec_name;
		    	size_setw_ = depth;

		    	if (it != type->begin_child_data()) {
			    	stream << ",";
				}
		    	doc_type_ = key_record;
		    	stream << endl;
		    	stream << setw((depth) * padding_size) << "";
		    	print(stream, &*it, depth);
		    	write_description(stream, it->description(), padding_size*size_setw_);
		    	doc_type_ = full_record;
		    	print(stream, &*it, depth);
		    }
		    //type->set_reference( refs[0] + refs[1] );
		    set_reference(type->type_name(), refs[0] + refs[1]);
		    remove_key(type->type_name());

		    //stream << endl;
			break;
	}

}


void OutputJSONTemplate::print_impl(ostream& stream, const Selection *type, unsigned int depth) {
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
		write_value(stream, value_);
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
        if (! type->made_extensive_doc()) {
            type->set_made_extensive_doc(true);

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
            stream << "{\\AddDoc{" << type->type_name() +"}}{"  << type->description() << "}";
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

                stream << "{\\AddDoc{" << type->type_name() << "::" << it->key_ << "}}{"
                       << it->description_ << "}" << endl;
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
        if (! type->made_extensive_doc()) {

            // Extensive description
            type->set_made_extensive_doc(true);

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
            stream << "{\\AddDoc{" << type->type_name() << "}}{"  << type->description() << "}" << endl;

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
        if (! type->made_extensive_doc()) {
            type->set_made_extensive_doc(true);

            stream <<endl << "\\begin{SelectionType}{" << internal::hyper_target("IT", type->type_name() ) << "}" <<endl;
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





std::ostream& operator<<(std::ostream& stream, OutputText type_output) {
    return type_output.print(stream) << endl;
}



std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output) {
    return type_output.print(stream) << endl;
}



std::ostream& operator<<(std::ostream& stream, OutputLatex type_output) {
    return type_output.print(stream) << endl;
}

} // closing namespace Type
} // closing namespace Input
