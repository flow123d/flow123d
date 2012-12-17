/*
 * type_output.hh
 *
 *  Created on: Nov 26, 2012
 *      Author: jb
 */

#ifndef TYPE_OUTPUT_HH_
#define TYPE_OUTPUT_HH_


#include "input/type_base.hh"
#include "input/type_record.hh"


namespace Input {

namespace Type {
/***
 * TODO:
 * - ve tride OutputBase navrhnout privatni virtualni medoty pro
 *
 */

class OutputBase {
public:
    OutputBase(TypeBase *type, unsigned int depth = 0) : type_(type), depth_(depth) {}

    void print(ostream& stream);

protected:
    /// Padding of new level of printout
    static const unsigned int padding_size = 4;

    /// Types of documentation output
    enum DocumentationType {
    	key_record,
    	full_record
    };

    // destructor
    virtual ~OutputBase();

    // data getters
    void get_array_sizes(Array array, unsigned int &lower , unsigned int &upper );
    void get_record_key(Record rec, unsigned int key_idx, Record::Key &key);
    void get_integer_bounds(Integer integer, int &lower , int &upper );
    void get_double_bounds(Double dbl, double &lower , double &upper );


    // type resolution like in json_to_storage
    void print(ostream& stream, const TypeBase *type);


    // following methods realize output in particular format
    // using getters from the base class OutputBase
    virtual void print(ostream& stream, const Record *type) = 0;
    virtual void print(ostream& stream, const Array *type) = 0;
    virtual void print(ostream& stream, const AbstractRecord *type) = 0;
    virtual void print(ostream& stream, const Selection *type) = 0;
	virtual void print(ostream& stream, const Integer *type) = 0;
	virtual void print(ostream& stream, const Double *type) = 0;
	virtual void print(ostream& stream, const Bool *type) = 0;
	virtual void print(ostream& stream, const String *type) = 0;
    virtual void print(ostream& stream, const FileName *type) = 0;

    void write_description(std::ostream& stream, const string& str);

    /// Object for which is created printout
    TypeBase *type_;
    /// Depth of printout
    unsigned int depth_;
    /// Actual level of printout
    unsigned int level_;
    /// Type of documentation output
    DocumentationType doc_type_;

};

class OutputText : public OutputBase {
public:
	OutputText(TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

	void print(ostream& stream) { OutputBase::print(stream); }

protected:

	void print(ostream& stream, const TypeBase *type) { OutputBase::print(stream, type); }

    void print(ostream& stream, const Record *type);
    void print(ostream& stream, const Array *type);
    void print(ostream& stream, const AbstractRecord *type);
    void print(ostream& stream, const Selection *type);
	void print(ostream& stream, const Integer *type);
	void print(ostream& stream, const Double *type);
	void print(ostream& stream, const Bool *type);
	void print(ostream& stream, const String *type);
    void print(ostream& stream, const FileName *type);
};

/**
 * Usage:
 * cout << OutputText( &my_record, 3) << "konec" << endl;
 */
std::ostream& operator<<(std::ostream& stream, OutputText type_output);



/*class OutputJSONTemplate : public OutputBase {

};


std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output);*/


} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_OUTPUT_HH_ */
