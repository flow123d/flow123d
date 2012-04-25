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
 * TODO: Better documentation

  Record         = MAP,        "key" = record R, vector V or scalar S
  Vector         = VECTOR,     "key" = vector[ (R or V or S), (R V S), (R V S), ... ]
  Value (scalar) = basic type  = string, number, bool, null

All classes: acquire subtree:
  val_type get_type()
  node &   get_key("key")
  node &   get_key("key", & def_tree )
  node &   get_item( index )
  node &   get_item( index, & def_tree )

Acquire scalar:
  get_string(), get_int(), get_double(), get_float(), get_bool()

  get_*( void ) - if not found, program fails

  get_*( & def_val ) - if not found, returns default

  get_*_check( & err_code ) - if not found, return non-zero error code and continue

Creating new nodes:
  Generic_node() (prev_node)
  Record_node() (prev_node)
  Vector_node() (prev_node)
  Value_node()(int)(double)(char*)(string)(bool) + variants s prev_node

Inserting values into tree:
  pro Record_node: void insert_key( "key", & node ) - insert; overwrite if exists
  pro Vector_node: void insert_item( index, & node ) - insert; overwrite if exists
  pro Value_node:  (int,double,char*,str,bool) set_value( ... ) - set value and return it
                   void  set_null()
Final usage:
  root.get_key("output").get_key("step").get_int();
  root.get_key("output").get_key("step").get_int(0);

  root.get_key("output").get_key("step").get_int_check( my_err );
  if ( my_err != EXIT_SUCCESS ) {
    //recover from fail here
  }

TODO:
* Pri vkladani insert_key(...) & insert_item(...) nastavovat prev_node
* mazani nodu & hodnot...
* get_parrent() (protected)
* presunout as_record(), as_value() a as_vector() do protected a zabudovat dovnitr get_*() metod - vzdy kontrolovat
  (nevadi, ze to bude pomalejsi, hlavni je, ze to bude jednodussi)
* value_type_to_string predelat na vector<string>, testovat delku vektoru
* gtest v test_units/input, make generic_node_test (generic_node_test.cpp)
* smazat include json.h z data_tree.h

Funkce data_tree:
    1) nacist JSON
    2) rozbalit REF
    3) s pomoci declaration_tree projit strom a vycistit ho - rovnou mazat nedeklarovane hodnoty

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
    static string  value_type_to_string[10];    //necessary for enum Value_type to string translation
                                                //lazy initialization = on first use
    static bool    value_type_to_string_filled; //already initialized?
protected:
    Value_type            value_type_;
    Generic_node &        prev_node_;

    //class data
    static Generic_node * empty_node_generic_; //empty instance
    static Record_node  * empty_node_record_;  //empty instance
    static Vector_node  * empty_node_vector_;  //empty instance
    static Value_node   * empty_node_value_;   //empty instance

    Generic_node():value_type_(type_generic),prev_node_(*this) {}
    Generic_node( Generic_node & prev_node ):value_type_(type_generic),prev_node_(prev_node) {}
    //Generic_node( Generic_node const & to_copy ); //copy constructor - implicit should be enough...
    Generic_node( const Value_type value_type ):value_type_(value_type),prev_node_( *this ) {};
    Generic_node( const Value_type value_type, Generic_node & prev_node ):value_type_(value_type),prev_node_(prev_node) {};

    virtual Generic_node & get_key( const string & key );
    virtual Generic_node & get_key( const string & key, Generic_node & default_tree );
    virtual Generic_node & get_key_check( const string & key, int & err_code );


public:
    /* Can not implement as_* here, need to know full class declaration.
     * Forward declaration is not enough (does not know available methods etc.)
     */
    virtual Record_node & as_record( void );
    virtual Vector_node & as_vector( void );
    virtual Value_node & as_value( void );

    Value_type get_type( void ) const { return value_type_; } //get node type
    const string & get_type_str( void ) const;                //get node type as a string description

    virtual Generic_node & get_item( const size_t id );
    virtual Generic_node & get_item( const size_t id, Generic_node & default_tree );
    virtual Generic_node & get_item_check( const size_t id, int & err_code );

    virtual size_t get_array_size( void );

    /*
     * Acquiring of final scalar values
     * w/o default value: fail if not found and no implicit conversion possible
     * w/  default value: always succeed, return default if not found or can not convert
     * w/  error code   : always succeed, return zero (or equivalent)
     *                    and set error code if not found or can not convert
     */

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

    virtual bool is_null( void );
    virtual bool not_null( void );

    // traverse recursively whole tree and print to stream
    friend ostream & operator<<( ostream & stream, Generic_node & node );

    virtual ~Generic_node();
};


} //namespace
#endif /* GENERIC_NODE_HPP_ */
