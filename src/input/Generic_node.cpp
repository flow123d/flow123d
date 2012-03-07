#include <string>
#include "Generic_node.hpp"
#include "Record_node.hpp"
#include "Vector_node.hpp"
#include "Value_node.hpp"

namespace flow {

//inicializace class dat
Record_node   * Generic_node::empty_node_record_  = new Record_node();  //prazdna instance
Vector_node   * Generic_node::empty_node_vector_  = new Vector_node();  //prazdna instance
Value_node    * Generic_node::empty_node_value_   = new Value_node();   //prazdna instance
Generic_node  * Generic_node::empty_node_generic_ = new Generic_node(); //prazdna instance
string          Generic_node::value_type_to_string[10];    //definice
bool            Generic_node::value_type_to_string_filled = false; //definice

const string & Generic_node::get_type_str( void ) const
{
    //lazy inicializace dat - az pri prvnim pouziti
    if ( !value_type_to_string_filled ) {
        value_type_to_string[type_generic] = "type_generic";
        value_type_to_string[type_string] = "type_string";
        value_type_to_string[type_number] = "type_number";
        value_type_to_string[type_record] = "type_record";
        value_type_to_string[type_vector] = "type_vector";
        value_type_to_string[type_bool] = "type_bool";
        value_type_to_string[type_null] = "type_null";
        value_type_to_string_filled = true;
    }

    return value_type_to_string[this->value_type_];
}

Generic_node & Generic_node::get_item(const size_t id)
{
    //v generic nic nemuze byt, vracime prazdnou
    return *empty_node_generic_;
}

Generic_node & Generic_node::get_key(const string & key)
{
    //v generic nic nemuze byt, vracime prazdnou
    return *empty_node_generic_;
}

Generic_node & Generic_node::get_item(const size_t id, Generic_node & default_tree)
{
    //v generic nic nemuze byt, vracime default
    return default_tree;
}

Generic_node & Generic_node::get_key(const string & key, Generic_node & default_tree)
{
    //v generic nic nemuze byt, vracime default
    return default_tree;
}

Record_node & Generic_node::as_record(void)
{
   if ( value_type_ == type_record ) {
       return * dynamic_cast < Record_node * > (this) ;
   } else {
       //zkousime pristoupit jako k recordu, ale neni record - vrat empty
       return *empty_node_record_;
   }
}

Vector_node & Generic_node::as_vector(void)
{
    if ( value_type_ == type_vector ) {
        return * dynamic_cast < Vector_node * > (this);
    } else {
        //zkousime pristoupit jako k vektoru, ale neni vektor - vrat empty
        return *empty_node_vector_;
    }
}

Value_node & Generic_node::as_value(void)
{
    if ( ( value_type_ == type_string ) || ( value_type_ == type_number ) ||
            ( value_type_ == type_bool ) || ( value_type_ == type_null ) ) {
        return * dynamic_cast < Value_node * > (this);
    } else {
        //zkousime pristoupit jako k value, ale neni value - vrat empty
        return *empty_node_value_;
    }
}

ostream & operator<<(ostream & stream, Generic_node & node)
{
    switch (node.get_type()) {
    case type_string:
    case type_number:
    case type_bool:
    case type_null:
        cout << node.as_value(); //spolecne pro vsechny skalarni hodnoty (dynamic cast)
        break;
    case type_record:
        cout << node.as_record(); //record vypis jako record (dynamic cast)
        break;
    case type_vector:
        cout << node.as_vector(); //vector vypis jako vector (dynamic cast)
        break;
    default:
        break;
    }
    return stream;
}

Generic_node & Generic_node::get_key_check(const string & key, int & err_code)
{
    //v generic nic nemuze byt, vracime prazdnou, s chybou
    err_code = 1;
    return *empty_node_generic_;
}

Generic_node & Generic_node::get_item_check(const size_t id, int & err_code)
{
    //v generic nic nemuze byt, vracime prazdnou, s chybou
    err_code = 1;
    return *empty_node_generic_;
}

bool Generic_node::get_bool(void) {
    //TODO: fail
    return false;
}

bool Generic_node::get_bool(const bool & default_value) {
    return default_value;
}

bool Generic_node::get_bool_check(int & err_code) {
    err_code = 1;
    return false;
}

int Generic_node::get_int(void) {
    //TODO: fail
    return 0;
}

int Generic_node::get_int(const int & default_value) {
    return default_value;
}

int Generic_node::get_int_check(int & err_code) {
    err_code = 1;
    return 0;
}

float Generic_node::get_float(void) {
    //TODO: fail
    return 0.0f;
}

float Generic_node::get_float(const float & default_value) {
    return default_value;
}

float Generic_node::get_float_check(int & err_code) {
    err_code = 1;
    return 0.0f;
}

double Generic_node::get_double(void) {
    //TODO: fail
    return 0.0;
}

double Generic_node::get_double(const double & default_value) {
    return default_value;
}

double Generic_node::get_double_check(int & err_code) {
    err_code = 1;
    return 0.0;
}

string Generic_node::get_string(void) {
    //TODO:fail
    return "";
}

string Generic_node::get_string(const string & default_value) {
    return default_value;
}

string Generic_node::get_string_check(int & err_code) {
    err_code = 1;
    return "";
}

Generic_node::~Generic_node() {
    //Tady nic.
    //Dynamicka data jsou spolecna pro tridu, nesmi se dealokovat.
}

}
