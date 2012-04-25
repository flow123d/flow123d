#include <string>
#include "../system/system.hh"
#include "Generic_node.hpp"
#include "Vector_node.hpp"
#include "Record_node.hpp"
#include "Value_node.hpp"

namespace flow {

//class data initialization
Record_node   * Generic_node::empty_node_record_  = new Record_node();  //empty instance
Vector_node   * Generic_node::empty_node_vector_  = new Vector_node();  //empty instance
Value_node    * Generic_node::empty_node_value_   = new Value_node();   //empty instance
Generic_node  * Generic_node::empty_node_generic_ = new Generic_node(); //empty instance
string          Generic_node::value_type_to_string[10];    //translation array definition
bool            Generic_node::value_type_to_string_filled = false; //lazy allocation helper definition

const string & Generic_node::get_type_str( void ) const
{
    //lazy data initialization - at the first use
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
    //generic is always empty - return empty
    return *empty_node_generic_;
}

Generic_node & Generic_node::get_key(const string & key)
{
    //generic is always empty - return empty
    return *empty_node_generic_;
}

Generic_node & Generic_node::get_item(const size_t id, Generic_node & default_tree)
{
    //generic is always empty - return default
    return default_tree;
}

Generic_node & Generic_node::get_key(const string & key, Generic_node & default_tree)
{
    //generic is always empty - return default
    return default_tree;
}

Record_node & Generic_node::as_record(void)
{
   if ( value_type_ == type_record ) {
       return * dynamic_cast < Record_node * > (this) ;
   } else {
       //wrong access as Record - return empty
       return *empty_node_record_;
   }
}

Vector_node & Generic_node::as_vector(void)
{
    if ( value_type_ == type_vector ) {
        return * dynamic_cast < Vector_node * > (this);
    } else {
        //wrong access as Vector - return empty
        return *empty_node_vector_;
    }
}

Value_node & Generic_node::as_value(void)
{
    if ( ( value_type_ == type_string ) || ( value_type_ == type_number ) ||
            ( value_type_ == type_bool ) || ( value_type_ == type_null ) ) {
        return * dynamic_cast < Value_node * > (this);
    } else {
        //wrong access as Value - return empty
        return *empty_node_value_;
    }
}

ostream & operator<<(ostream & stream, Generic_node & node)
{
    switch (node.value_type_) {
    case type_string:
    case type_number:
    case type_bool:
    case type_null:
        stream << node.as_value(); //common for all scalar types (dynamic cast)
        break;
    case type_record:
        stream << node.as_record(); //Record print as Record (dynamic cast)
        break;
    case type_vector:
        stream << node.as_vector(); //Vector print as Vector (dynamic cast)
        break;
    default:
        stream << "Err: instance of Generic_node should not exist.";
        break;
    }
    return stream;
}

Generic_node & Generic_node::get_key_check(const string & key, int & err_code)
{
    //generic is always empty - return empty and error
    err_code = 1;
    return *empty_node_generic_;
}

Generic_node & Generic_node::get_item_check(const size_t id, int & err_code)
{
    //generic is always empty - return empty and error
    err_code = 1;
    return *empty_node_generic_;
}

bool Generic_node::get_bool(void) {
    xprintf(PrgErr, "Can not get_bool() from Generic_node." );
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
    xprintf(PrgErr, "Can not get_int() from Generic_node." );
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
    xprintf(PrgErr, "Can not get_float() from Generic_node." );
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
    xprintf(PrgErr, "Can not get_double() from Generic_node." );
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
    xprintf(PrgErr, "Can not get_string() from Generic_node." );
    return "";
}

string Generic_node::get_string(const string & default_value) {
    return default_value;
}

string Generic_node::get_string_check(int & err_code) {
    err_code = 1;
    return "";
}

void Generic_node::delete_id(const size_t id) {
    xprintf(PrgErr, "Can not delete_id(), not in Vector_node." );
}

void Generic_node::delete_key(const string& key) {
    xprintf(PrgErr, "Can not delete_key(), not in Record_node." );
}

size_t Generic_node::get_array_size(void) {
    return 0;
}

bool Generic_node::is_null(void) {
    return true;
}

bool Generic_node::not_null(void) {
    return false;
}

Generic_node::~Generic_node() {
    //The only dynamic data present in Generic_node are Class data - must not be deallocated.
}

} //namespace

