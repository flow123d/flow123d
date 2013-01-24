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
 * - v OutputText: Record ktery dedi bude uvadet jmeno predka (AbstractRecordu) - TODO: pridat shared_ptr do Record ukazujici na predka (trvale)
 *
 * - nejak informovat o povoleni automaticke konverze pro Record
 *
 * - OutputJSONTemplate:
 *   - vypis defaultnich hodnot pro skalarni klice, typ klice do komentare na stejnem radku
 *     pokud je obligatory vypsat <OBLIGATORY>, klic nezakomentovat
 *     Optional: prefix klice "OPT_", nezakomentovat (funguje jako komentar.. neznamy klic), jako hodnotu pouzit Integer, Double:0, String, Selection: ""
 *     Read_time: prefix "OPT_" , za rovnitko komentar?? co s default at read time a optional
 *     patrne je uvest ale v komentari, idea je aby soubor byl platny JSON soubor az na OBLIGATORY klice, ktere uzivatel musi vyplnit
 *
 *   - Pro selection uvest mozne hodnoty a jejich popisy do komentare
 *     # Selection of 3 values:
       # "None" - Mortar space: P0 on elements of lower dimension.
       # "P0"   - Mortar space: P0 on elements of lower dimension.
       # "P1"   - Mortar space: P1 on intersections, using non-conforming pressures.
       # ---------
       # Method for coupling Darcy flow between dimensions.
       mortar_method = "None"
 *
 *   - Pro abstract record: uvest TYPE (selection), pak jeho klice (dedi se) a pak
 *     jednotlive potomky pomoci klicu s pripojenym typem recordu,
 *   - pro klice ktere maji typ ktery uz byl popsan se uvede reference na klic kde se popis vyskytuje
 *
 *     # abstract record FieldBase_3_to_1x1_double
 *     # description:
 *     # ...
 *     # ----------------------------------------------- DESCENDANTS FOLLOWS
 *     # record FieldConstant, descendant of FieldBase_3_to_1x1_double
 *     init_pressure={
 *          TYPE="FieldConstant"
 *
 *          # description
 *          value=
 *     },
 *     # record FieldFormula
 *     init_pressure={
 *          TYPE="FieldFormula"
 *          formula=
 *          parameters=
 *     }
 *
 *     # abstract record FieldBase_3_to_1x1_double
 *     # description:
 *     water_source={REF="/.../init_pressure"},
 *
 *     # Array, size limits [2,3]
 *     # key description ...
 *     array_key=[
 *         # Record ...
 *         #
 *         {
 *             ...
 *         },
 *         < 1 more entry >
 *         ]
 *
 *
 *   - ?? jak se vyporadat s default values  pro automaticky konvertovatelne recordy a pole?
 *   - avoid repetitive output of same tyes
 *
 *   - jak je reseno potlaceni opakovaneho vypisovani Recordu, je mozno odstranit reset_doc_flags z TypeBase
 *   - pri dalsich upravach v exceptions zkusit predavat v EI smart_ptr na TypeBase, nicmene na miste vyvolani fce by muselo byt
 *     make_shared, aby se vytvorila kopie spravneho typu. (faktycky je potreba objekt Type kopirovat jinak nevim, ze mi ho nekdo neodalokuje.
 */

class OutputBase {
public:
    OutputBase(const TypeBase *type, unsigned int depth = 0);

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

    /**
     * Write out a string with given padding of every new line.
     *
     * @param stream Output stream
     * @param str Printed description
     * @param hash_count Count of '#' chars in description
     */
    virtual void write_description(std::ostream& stream, const string& str, unsigned int hash_count = 1) = 0;

    void write_value(std::ostream& stream, Default dft);

    /// Object for which is created printout
    const TypeBase *type_;
    /// Depth of printout
    unsigned int depth_;
    /// Type of documentation output
    DocumentationType doc_type_;

    /// temporary value for printout of description (used in std::setw function)
    unsigned int size_setw_;

};

class OutputText : public OutputBase {
public:
	OutputText(const TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

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

    void write_description(std::ostream& stream, const string& str, unsigned int hash_count = 1);

};

/**
 * Usage:
 * cout << OutputText( &my_record, 3) << "konec" << endl;
 */
std::ostream& operator<<(std::ostream& stream, OutputText type_output);



class OutputJSONTemplate : public OutputBase {
public:
	OutputJSONTemplate(TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

	void print(ostream& stream) {
		key_name_ = "";
		OutputBase::print(stream);
	}

protected:
	void print(ostream& stream, const TypeBase *type, unsigned int depth = 0) {
		OutputBase::print(stream, type, depth);
	}

    void print(ostream& stream, const Record *type, unsigned int depth = 0);
    void print(ostream& stream, const Array *type, unsigned int depth = 0);
    void print(ostream& stream, const AbstractRecord *type, unsigned int depth = 0);
    void print(ostream& stream, const Selection *type, unsigned int depth = 0);
	void print(ostream& stream, const Integer *type, unsigned int depth = 0);
	void print(ostream& stream, const Double *type, unsigned int depth = 0);
	void print(ostream& stream, const Bool *type, unsigned int depth = 0);
	void print(ostream& stream, const String *type, unsigned int depth = 0);
    void print(ostream& stream, const FileName *type, unsigned int depth = 0);

    void write_description(std::ostream& stream, const string& str, unsigned int hash_count = 1);

private:
    /**
     * Prints value according to DefaultType
     * Respects obligatory, optional and read time flag
     *
     * @param stream Output stream
     * @param depth Depth of output
     * @param empty_val Default empty value (zero for numeric types, empty string ...)
     * @param invalid_val Flag if value is invalid for its type
     * @param has_quote Flag if value is enclosed in quotes
     */
    void print_default_value(ostream& stream, unsigned int depth, string empty_val, bool invalid_val, bool has_quote = false);

    /// temporary value of actually record type
    string key_name_;
    /// temporary value of actually record description
    string description_;
    /// temporary value of actually record value
    Default value_;
};


std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output);


} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_OUTPUT_HH_ */
