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


void OutputBase::get_array_sizes(Array array, int &lower , int &upper ) {
	lower = array.lower_bound_;
	upper = array.upper_bound_;
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


void OutputBase::print(ostream& stream, const TypeBase *type) {
	if (typeid(*type) == typeid(Type::Record)) {
		print(stream, static_cast<const Type::Record *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Array)) {
		print(stream, static_cast<const Type::Array *>(type) );
	} else
	if (typeid(*type) == typeid(Type::AbstractRecord)) {
		print(stream, static_cast<const Type::AbstractRecord *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Selection)) {
		print(stream, static_cast<const Type::Selection *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Integer)) {
		print(stream, static_cast<const Type::Integer *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Double)) {
		print(stream, static_cast<const Type::Double *>(type) );
	} else
	if (typeid(*type) == typeid(Type::Bool)) {
		print(stream, static_cast<const Type::Bool *>(type) );
	} else {
		const Type::String * string_type = dynamic_cast<const Type::String *>(type);
        if (string_type != NULL ) {
        	print(stream, string_type );
        	return;
        }

		const Type::FileName * file_name_type = dynamic_cast<const Type::FileName *>(type);
        if (file_name_type != NULL ) {
        	print(stream, file_name_type );
        	return;
        }

        // default -> error
        xprintf(Err,"Unknown descendant of TypeBase class, name: %s\n", typeid(type).name());
	}
}




/*******************************************************************
 * implementation of OutputText
 */

void OutputText::print(ostream& stream, const Record *type) {
	// Print documentation of record
}


void OutputText::print(ostream& stream, const Array *type) {
	int lower_size, upper_size;

	get_array_sizes(*type, lower_size, upper_size);
	stream << "Array, size limits: [" << lower_size << ", " << upper_size << "] of type: " << endl;
	stream << setw(4) << "";
	//print(stream, type->get_sub_type());
}


void OutputText::print(ostream& stream, const AbstractRecord *type) {
	// Print documentation of abstract record
}


void OutputText::print(ostream& stream, const Selection *type) {
    /*stream << endl << "Selection '" << type->type_name() << "' of " << type->size() << " values." << endl;
    stream << setw(0) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
    // keys
    for (Selection::SelectionData::keys_const_iterator it = type->data_.keys_.begin(); it != type->data_.keys_.end(); ++it) {
        //stream << setw(4) << "" << it->key_ << " = " << it->value;
        //if (it->description_ != "")
        //    stream << " (" << it->description_ << ")";
        //stream << endl;
    }
    stream << setw(0) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type_name_ << endl;*/
}


void OutputText::print(ostream& stream, const Integer *type) {
	int lower_bound, upper_bound;
	get_integer_bounds(*type, lower_bound, upper_bound);
	stream << "Integer in [" << lower_bound << ", " << upper_bound << "]";
}


void OutputText::print(ostream& stream, const Double *type) {
	double lower_bound, upper_bound;
	get_double_bounds(*type, lower_bound, upper_bound);
	stream << "Double in [" << lower_bound << ", " << upper_bound << "]";
}


void OutputText::print(ostream& stream, const Bool *type) {
	stream << "Bool";
}


void OutputText::print(ostream& stream, const String *type) {
	stream << "String (generic)";
}


void OutputText::print(ostream& stream, const FileName *type) {
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


std::ostream& operator<<(std::ostream& stream, OutputText type_output) {
	//cout << type_output.print(stream) << endl;
	return stream;
}



} // closing namespace Type
} // closing namespace Input
