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
 * - proc je v OutputText pretizena resolution metoda print ?
 *   podobne v OutputJSONTemplate
 * - v OutputText: Record ktery dedi bude uvadet jmeno predka (AbstractRecordu)
 * - nejak informovat o povoleni automaticke konverze pro Array a Record
 * - OutptuJSONTemplate:
 *   - vypis defaultnich hodnot pro skalarni klice, typ klice do komentare na stejnem radku
 *     pokud je obligatory vypsat <OBLIGATORY>, ?? co s default at read time a optional
 *     patrne je uvest ale v komentari, idea je aby soubor byl platny JSON soubor az na OBLIGATORY klice, ktere uzivatel musi vyplnit
 *   - Pro selection uvest mozne hodnoty a jejich popisy do komentare
 *   - Pro abstract record: uvest TYPE (selection), pak jeho klice (dedi se) a pak
 *     jednotlive potomky pomoci klicu s pripojenym typem recordu,
 *   - pro klice ktere maji typ ktery uz byl popsan se uvede reference na klic kde se popis vyskytuje
 *
 *     # abstract record FieldBase_3_to_1x1_double
 *     # description:
 *     # ...
 *     init_pressure={
 *         # selection of 5 values
 *         # FieldConstant - ...
 *         # FieldFormula - ...
 *         # FieldInterpolationP0
 *         # FieldElementwise
 *         # FieldPython
 *         TYPE="FieldConstant"
 *
 *         # .. some common key, with default value 0
 *         common_key=0
 *     }
 *     # record FieldConstant, descendant of FieldBase_3_to_1x1_double
 *     init_pressure_FieldConstant={
 *     ...
 *     }
 *     # record FieldFormula
 *     init_pressure_FieldFormula={
 *     ...
 *     }
 *
 *     water_source={REF=/.../init_pressure}
 *   - ?? jak se vyporadat s default values  pro automaticky konvertovatelne recordy a pole?
 *   - avoid repetitive output of same tyes
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
    void print(ostream& stream, const TypeBase *type, unsigned int depth = 0);


    // following methods realize output in particular format
    // using getters from the base class OutputBase
    virtual void print(ostream& stream, const Record *type, unsigned int depth = 0) = 0;
    virtual void print(ostream& stream, const Array *type, unsigned int depth = 0) = 0;
    virtual void print(ostream& stream, const AbstractRecord *type, unsigned int depth = 0) = 0;
    virtual void print(ostream& stream, const Selection *type, unsigned int depth = 0) = 0;
	virtual void print(ostream& stream, const Integer *type, unsigned int depth = 0) = 0;
	virtual void print(ostream& stream, const Double *type, unsigned int depth = 0) = 0;
	virtual void print(ostream& stream, const Bool *type, unsigned int depth = 0) = 0;
	virtual void print(ostream& stream, const String *type, unsigned int depth = 0) = 0;
    virtual void print(ostream& stream, const FileName *type, unsigned int depth = 0) = 0;

    void write_description(std::ostream& stream, const string& str);

    /// Object for which is created printout
    TypeBase *type_;
    /// Depth of printout
    unsigned int depth_;
    /// Type of documentation output
    DocumentationType doc_type_;

};

class OutputText : public OutputBase {
public:
	OutputText(TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

	void print(ostream& stream) { OutputBase::print(stream); }

protected:

	void print(ostream& stream, const TypeBase *type, unsigned int depth = 0) { OutputBase::print(stream, type, depth); }

    void print(ostream& stream, const Record *type, unsigned int depth = 0);
    void print(ostream& stream, const Array *type, unsigned int depth = 0);
    void print(ostream& stream, const AbstractRecord *type, unsigned int depth = 0);
    void print(ostream& stream, const Selection *type, unsigned int depth = 0);
	void print(ostream& stream, const Integer *type, unsigned int depth = 0);
	void print(ostream& stream, const Double *type, unsigned int depth = 0);
	void print(ostream& stream, const Bool *type, unsigned int depth = 0);
	void print(ostream& stream, const String *type, unsigned int depth = 0);
    void print(ostream& stream, const FileName *type, unsigned int depth = 0);
};

/**
 * Usage:
 * cout << OutputText( &my_record, 3) << "konec" << endl;
 */
std::ostream& operator<<(std::ostream& stream, OutputText type_output);



class OutputJSONTemplate : public OutputBase {
public:
	OutputJSONTemplate(TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

	void print(ostream& stream) { OutputBase::print(stream); }

protected:
	void print(ostream& stream, const TypeBase *type, unsigned int depth = 0) { OutputBase::print(stream, type, depth); }

    void print(ostream& stream, const Record *type, unsigned int depth = 0);
    void print(ostream& stream, const Array *type, unsigned int depth = 0);
    void print(ostream& stream, const AbstractRecord *type, unsigned int depth = 0);
    void print(ostream& stream, const Selection *type, unsigned int depth = 0);
	void print(ostream& stream, const Integer *type, unsigned int depth = 0);
	void print(ostream& stream, const Double *type, unsigned int depth = 0);
	void print(ostream& stream, const Bool *type, unsigned int depth = 0);
	void print(ostream& stream, const String *type, unsigned int depth = 0);
    void print(ostream& stream, const FileName *type, unsigned int depth = 0);

};


std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output);


} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_OUTPUT_HH_ */
