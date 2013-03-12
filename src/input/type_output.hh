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
 *
 *
 *
 *  value = <OBLIGATORY>
 *          #    is String (generic)
            # String, array of strings, or matrix of strings with formulas for individual entries of scalar, vector, or tensor value respectively.
            # For vector values, you can use just one string to enter homogeneous vector.
            # For square NxN-matrix values, you can use:
            # * array of strings of size N to enter diagonal matrix
            # * array of strings of size (N+1)*N/2 to enter symmetric matrix (upper triangle, row by row)
            # * just one string to enter (spatially variable) multiple of the unit matrix.
            # Formula can contain variables x,y,z,t and usual operators and functions.
 *
 */

/**
 * Base abstract class for output description of the Input::Type tree.
 * Output into various formats is implemented by derived classes.
 *
 * Usage:
 * cout << OutputText( &my_record, 3) << "konec" << endl;
 *
 */
class OutputBase {
public:

    /**
     * Performs output of the documentation into given @p stream.
     * Returns reference to the same stream.
     */
    virtual ostream& print(ostream& stream);

protected:
    /// Types of documentation output
    enum DocumentationType {
        key_record,
        full_record
    };

    /// Padding of new level of printout, used where we use indentation.
    static const unsigned int padding_size = 4;


    /**
     * Constructor
     *
     * @param type Stores input sequence
     * @param depth Depth of output
     */
    OutputBase(const TypeBase *type, unsigned int depth = 0);


    // destructor
    virtual ~OutputBase();

    // data getters
    void get_array_sizes(Array array, unsigned int &lower , unsigned int &upper );
    void get_record_key(Record rec, unsigned int key_idx, Record::Key &key);
    void get_integer_bounds(Integer integer, int &lower , int &upper );
    void get_double_bounds(Double dbl, double &lower , double &upper );


    /**
     * Perform resolution according to actual @p type (using typeid) and call particular print_impl method.
     */
    void print(ostream& stream, const TypeBase *type, unsigned int depth);


    /**
     *  following methods realize output in particular format
     *  using getters from the base class OutputBase
     */
    virtual void print_impl(ostream& stream, const Record *type, unsigned int depth) = 0;
    virtual void print_impl(ostream& stream, const Array *type, unsigned int depth) = 0;
    virtual void print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth) = 0;
    virtual void print_impl(ostream& stream, const Selection *type, unsigned int depth) = 0;
	virtual void print_impl(ostream& stream, const Integer *type, unsigned int depth) = 0;
	virtual void print_impl(ostream& stream, const Double *type, unsigned int depth) = 0;
	virtual void print_impl(ostream& stream, const Bool *type, unsigned int depth) = 0;
	virtual void print_impl(ostream& stream, const String *type, unsigned int depth) = 0;
    virtual void print_impl(ostream& stream, const FileName *type, unsigned int depth) = 0;

    /**
     * Write out a string with given padding of every new line.
     *
     * @param stream Output stream
     * @param str Printed description
     * @param hash_count Count of '#' chars in description
     */
    //virtual void write_description(std::ostream& stream, const string& str, unsigned int hash_count = 1) = 0;
    /**
     * Output indented multi-line string.
     */
    void write_description(std::ostream& stream, const string& str, unsigned int padding, unsigned int hash_count = 1);


    /**
     * Write value stored in dft.
     *
     * Enclose value in quotes if it's needed or write info that value is optional or obligatory.
     */
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


/**********************************************************************************************************************/

/**
 * Class for create text documentation
 */
class OutputText : public OutputBase {
public:
	OutputText(const TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

protected:

    void print_impl(ostream& stream, const Record *type, unsigned int depth);
    void print_impl(ostream& stream, const Array *type, unsigned int depth);
    void print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const Selection *type, unsigned int depth);
	void print_impl(ostream& stream, const Integer *type, unsigned int depth);
	void print_impl(ostream& stream, const Double *type, unsigned int depth);
	void print_impl(ostream& stream, const Bool *type, unsigned int depth);
	void print_impl(ostream& stream, const String *type, unsigned int depth);
    void print_impl(ostream& stream, const FileName *type, unsigned int depth);


};







/**
 * Class for create and JSON template documentation
 */
class OutputJSONTemplate : public OutputBase {
public:
    /**
     * Constructor for output of the input type tree with root @p type.
     * The input type tree is searched by DFS algorithm into @p depth.
     */
	OutputJSONTemplate(TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

	/**
	 * Perform output of the documentation into given stream.
	 */
	ostream& print(ostream& stream);

protected:
	// Need to implement the resolution function. Just call that in the base class.
	void print(ostream& stream, const TypeBase *type, unsigned int depth) {
		OutputBase::print(stream, type, depth);
	}


    void print_impl(ostream& stream, const Record *type, unsigned int depth);
    void print_impl(ostream& stream, const Array *type, unsigned int depth);
    void print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const Selection *type, unsigned int depth);
	void print_impl(ostream& stream, const Integer *type, unsigned int depth);
	void print_impl(ostream& stream, const Double *type, unsigned int depth);
	void print_impl(ostream& stream, const Bool *type, unsigned int depth);
	void print_impl(ostream& stream, const String *type, unsigned int depth);
    void print_impl(ostream& stream, const FileName *type, unsigned int depth);

    //void write_description(std::ostream& stream, const string& str, unsigned int hash_count = 1);

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
    /// temporary value of actually record value
    Default value_;
};




/**
 * Class for create and Latex documentation
 */
class OutputLatex : public OutputBase {
public:
    OutputLatex(TypeBase *type, unsigned int depth = 0) : OutputBase(type, depth) {}

    ostream & print(ostream& stream);

protected:
    // Need to implement the resolution function. Just call that in the base class.
    void print(ostream& stream, const TypeBase *type, unsigned int depth) {
        OutputBase::print(stream, type, depth);
    }

    void print_impl(ostream& stream, const Record *type, unsigned int depth);
    void print_impl(ostream& stream, const Array *type, unsigned int depth);
    void print_impl(ostream& stream, const AbstractRecord *type, unsigned int depth);
    void print_impl(ostream& stream, const Selection *type, unsigned int depth);
    void print_impl(ostream& stream, const Integer *type, unsigned int depth);
    void print_impl(ostream& stream, const Double *type, unsigned int depth);
    void print_impl(ostream& stream, const Bool *type, unsigned int depth);
    void print_impl(ostream& stream, const String *type, unsigned int depth);
    void print_impl(ostream& stream, const FileName *type, unsigned int depth);


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
    //string key_name_;
    /// temporary value of actually record value
    //Default value_;
};


/**
 * Overrides output operator for simple output of the input type tree.
 */
std::ostream& operator<<(std::ostream& stream, OutputText type_output);
std::ostream& operator<<(std::ostream& stream, OutputJSONTemplate type_output);
std::ostream& operator<<(std::ostream& stream, OutputLatex type_output);



} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_OUTPUT_HH_ */
