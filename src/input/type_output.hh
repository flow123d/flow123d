/*
 * type_output.hh
 *
 *  Created on: Nov 26, 2012
 *      Author: jb
 */

#ifndef TYPE_OUTPUT_HH_
#define TYPE_OUTPUT_HH_


namespace Input {

namespace Type {
/***
 * TODO:
 * - ve tride OutputBase navrhnout privatni virtualni medoty pro
 *
 */

class OutputBase {
public:
    OutBase(TypeBase *type, unsigned int depth = 0);
private:
    // data getters
    void get_array_sizes(Array array, int &lower , int  &upper );
    void get_record_key(Record rec, unsigned int key_idx, Record::Key &key)


    // type resolution like in json_to_storage
    void print(ostream& stream, const TypeBase *type);


    // following methods realize output in particular format
    // using getters from the base class OutputBase
    virtual void print(ostream& stream, const Record *type);
    virtual void print(ostream& stream, const Array *type);
    virtual void print(ostream& stream, const AbstractRecord *type);
    virtual void print(ostream& stream, const Selection *type);
    virtual void print(ostream& stream, const Scalar *type);
    virtual void print(ostream& stream, const FileName *type);

    unsigned int detph_

};

class OutputText : public OutputBase {

};

/**
 * Usage:
 * cout << OutputText( &my_record, 3) << "konec" << endl;
 */
std::ostream& operator<<(std::ostream& stream, OutputText type_output);



class OutputJSONTemplate : public OutputBase {

};


std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output);


} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_OUTPUT_HH_ */
