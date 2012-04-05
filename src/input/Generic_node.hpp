#ifndef GENERIC_NODE_HPP_
#define GENERIC_NODE_HPP_

#include <map>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

namespace flow {

class Record_node;
class Vector_node;
class Value_node;

/*
 * TODO: PRELOZIT, DOKUMENTACE
Kam to smeruje:

  Record         = cca MAP,        "klic"= record R , vektor V nebo skalar S
  Vector         = cca VECTOR,     "klic"= vector[ (R nebo V nebo S), (R V S), (R V S), ... ]
  Value (skalar) = cca basic type  = string, number, bool, null

Pro vsechny tridy, ziskani podstromu:
  val_type get_type()
  node &   get_key("klic")
  node &   get_key("klic", & def_tree )
  node &   get_item( index )
  node &   get_item( index, & def_tree )

Moznosti ziskani skalarnich hodnot:
  get_string(), get_int(), get_double(), get_float(), get_bool()

  get_*( void ) - kdyz neni, pada

  get_*( & def_val ) - kdyz neni, vrati default

  get_*_check( & err_code ) - kdyz neni, vrati chybu a pokracuje

Vytvareni novych nodu:
  Generic_node() (prev_node)
  Record_node() (prev_node)
  Vector_node() (prev_node)
  Value_node()(int)(double)(char*)(string)(bool) + varianty s prev_node

Moznosti vkladani hodnot:
  pro Record_node: void insert_key( "klic", & node ) - vlozi; pokud existuje, prepise
  pro Vector_node: void insert_item( index, & node ) - vlozi; pokud existuje, prepise
  pro Value_node:  (int,double,char*,str,bool) set_value( ... ) - nastavi hodnotu a rovnou ji vrati
                   void  set_null()
Finalni pouziti:
  root.get_key("output").get_key("step").get_int();
  root.get_key("output").get_key("step").get_int(0);

  root.get_key("output").get_key("step").get_int_check( my_err );
  if ( my_err != EXIT_SUCCESS ) {
    //zde osetrim chybku
  }

TODO:
* Pri vkladani insert_key(...) & insert_item(...) nastavovat prev_node

* Pri chybe:
* anonymous_node (potomek generic_node)
pamatuje posledni validni generic node ve stromu a ma zasobnik operaci co se po nem chtelo
umi :
1) vytvorit validni node vcetne nadrazenych nodu pri deklaraci:

   anon_node=parent.get_key("key_not_on_input");
   anon_node.get_key("some_key_1").declare(Bool, false, "Description");
   anon_node.get_key("some_key_2").declare(Int, 1, "Description");

2) vypisy chyb, pri cteni value bez default hodnoty

*/


/*!
 * Types of nodes in graph
 */
enum Value_type { type_generic, type_record, type_vector, type_string, type_number, type_bool, type_null };

string & Value_type_to_str( const Value_type vt );

/*!
 * @brief Generic node - interface class for JSON input layer
 *
 * This class specify interface for reading values through the JSON input layer. The input layer consists
 * from tree of input nodes of different types but all of them are derived form Generic_node since need to store some of them in common container.
 * This implies that Generic_node provides all methods
 *
 * The leaves of this tree are instances of derived class Value_node that represents any
 * basic JSON type namely string, number, boolean or NULL.  The class Value_node overrides methodsWe provide chich is the root class
 */
class Generic_node {
private:
    static string  value_type_to_string[10];    //pro preklad enum Value_type na string.
                                                //Lazy inicializace = pri prvnim cteni.
    static bool    value_type_to_string_filled; //Uz je inicializovano?
protected:
    Value_type            value_type_;
    Generic_node &        prev_node_;

    //class data
    static Generic_node * empty_node_generic_; //prazdna instance
    static Record_node  * empty_node_record_;  //prazdna instance
    static Vector_node  * empty_node_vector_;  //prazdna instance
    static Value_node   * empty_node_value_;   //prazdna instance

    //constructor s urcenim datoveho typu - pristupny pouze z potomku
    Generic_node( const Value_type value_type ):value_type_(value_type),prev_node_( *this ) {};
    Generic_node( const Value_type value_type, Generic_node & prev_node ):value_type_(value_type),prev_node_(prev_node) {};

public:
    Generic_node():value_type_(type_generic),prev_node_(*this) {}
    Generic_node( Generic_node & prev_node ):value_type_(type_generic),prev_node_(prev_node) {}
    //Generic_node( Generic_node const & to_copy ); //copy constructor - staci implicitni...

    Value_type get_type( void ) const { return value_type_; } //ziska typ nodu. To by mel umet kazdy
    const string & get_type_str( void ) const;                //ziska typ nodu jako string popis

    virtual Generic_node & get_key( const string & key );
    virtual Generic_node & get_key( const string & key, Generic_node & default_tree );
    virtual Generic_node & get_key_check( const string & key, int & err_code );
    virtual Generic_node & get_item( const size_t id );
    virtual Generic_node & get_item( const size_t id, Generic_node & default_tree );
    virtual Generic_node & get_item_check( const size_t id, int & err_code );

    /* Implementace as_* tady nefunguje, je potreba znat presnou deklaraci tridy.
     * Forward deklarace nestaci (nevi, jake ma k dispozici metody apod.)
     */
    virtual Record_node & as_record( void );
    virtual Vector_node & as_vector( void );
    virtual Value_node & as_value( void );

    //ziskani uz finalnich skalarnich hodnot
    //bez def: pri nemoznosti konverze spadne
    //  s def: uspeje vzdy, protoze pri nemoznosti konverze pouzije default
    //  s err: uspeje vzdy, pri neuspechu vyplni err, jako hodnotu vrati 0, null, nebo ekvivalent
    virtual bool get_bool( void );
    virtual bool get_bool( const bool & default_value );
    virtual bool get_bool_check( int & err_code );

    virtual int get_int( void );
    virtual int get_int( const int & default_value );
    virtual int get_int_check( int & err_code );

    virtual float get_float( void );
    virtual float get_float( const float & default_value );
    virtual float get_float_check( int & err_code );

    virtual double get_double( void );
    virtual double get_double( const double & default_value );
    virtual double get_double_check( int & err_code );

    virtual string get_string( void );
    virtual string get_string( const string & default_value );
    virtual string get_string_check( int & err_code );


    // hierarchicky projde cely podstrom do hloubky a vypise ho
    friend ostream & operator<<( ostream & stream, Generic_node & node );

    virtual ~Generic_node();
};


} //namespace
#endif /* GENERIC_NODE_HPP_ */
