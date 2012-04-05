#include "Value_node.hpp"

namespace flow {

ostream & operator <<(ostream & stream, Value_node & node)
{
    switch ( node.value_type_ ) {
    case type_string:
        stream << "\"" << node.value_string_ << "\"";
        break;
    case type_bool:
        stream << (node.value_bool_?"true":"false");
        break;
    case type_number:
        stream << node.value_number_;
        break;
    case type_null:
        stream << "null";
        break;
    default:
        break;
    }

    return stream;
}

bool Value_node::get_bool(void) {
    if ( value_type_ == type_bool )
        return value_bool_;
    else
        return Generic_node::get_bool();
}

bool Value_node::get_bool(const bool & default_value) {
    if ( value_type_ == type_bool )
        return value_bool_;
    else
        return default_value;
}

bool Value_node::get_bool_check(int & err_code) {
    if ( value_type_ == type_bool )
    {
        err_code = 0;
        return value_bool_;
    }
    else
        return Generic_node::get_bool_check(err_code);
}

int Value_node::get_int(void) {
    if ( value_type_ == type_number )
        return (int) (value_number_);
    else
        return Generic_node::get_int();
}

int Value_node::get_int(const int & default_value) {
    if ( value_type_ == type_number )
        return (int) (value_number_);
    else
        return default_value;
}

int Value_node::get_int_check(int & err_code) {
    if ( value_type_ == type_number )
    {
        err_code = 0;
        return (int) (value_number_);
    }
    else
        return Generic_node::get_int_check(err_code);
}

float Value_node::get_float(void) {
    if ( value_type_ == type_number )
        return (float)(value_number_);
    else
        return Generic_node::get_float();
}

float Value_node::get_float(const float & default_value) {
    if ( value_type_ == type_number )
        return (float)(value_number_);
    else
        return default_value;
}

float Value_node::get_float_check(int & err_code) {
    if ( value_type_ == type_number )
    {
        err_code = 0;
        return (float)(value_number_);
    }
    else
        return Generic_node::get_float_check(err_code);
}

double Value_node::get_double(void) {
    if ( value_type_ == type_number )
        return (double)(value_number_);
    else
        return Generic_node::get_double();
}

double Value_node::get_double(const double & default_value) {
    if ( value_type_ == type_number )
        return (double)(value_number_);
    else
        return default_value;
}

double Value_node::get_double_check(int & err_code) {
    if ( value_type_ == type_number )
    {
        err_code = 0;
        return (double)(value_number_);
    }
    else
        return Generic_node::get_double_check(err_code);
}

string Value_node::get_string(void) {
    if ( value_type_ == type_string )
        return value_string_;
    else
        return Generic_node::get_string();
}

string Value_node::get_string(const string & default_value) {
    if ( value_type_ == type_string )
        return value_string_;
    else
        return default_value;
}

string Value_node::get_string_check(int & err_code) {
    if ( value_type_ == type_string )
    {
        err_code = 0;
        return value_string_;
    }
    else
        return Generic_node::get_string_check(err_code);
}

Value_node::~Value_node() {
    //Empty.
}

}
