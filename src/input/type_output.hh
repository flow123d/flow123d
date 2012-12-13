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
    // destructor
    virtual ~OutputBase();

    // data getters
    void get_array_sizes(Array array, int &lower , int &upper );
    void get_record_key(Record rec, unsigned int key_idx, Record::Key &key);
    void get_integer_bounds(Integer integer, int &lower , int &upper );
    void get_double_bounds(Double dbl, double &lower , double &upper );
    //void get_selection_keys(Selection sel, double &lower , double &upper );


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

    TypeBase *type_;
    unsigned int depth_;

};

class OutputText : public OutputBase {
public:
	OutputText(TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

	void print(ostream& stream) { OutputBase::print(stream); }

protected:

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
