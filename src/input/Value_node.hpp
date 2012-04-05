#ifndef VALUE_NODE_HPP_
#define VALUE_NODE_HPP_

#include "Generic_node.hpp"

namespace flow {

/*!
 * @brief Scalar node - represents JSON construct "value",
 *        that can be "number, bool, string, null"
 */
class Value_node : public Generic_node {
    string              value_string_;
    double              value_number_;
    bool                value_bool_;

public:
    Value_node():Generic_node(type_null)                 {}
    Value_node( bool b ):Generic_node(type_bool)         { value_bool_ = b;}
    Value_node( int i ):Generic_node(type_number)        { value_number_ = i; }
    Value_node( double lf ):Generic_node(type_number)    { value_number_ = lf; }
    Value_node( char * str ):Generic_node(type_string)   { value_string_ = str; }
    Value_node( string & str ):Generic_node(type_string) { value_string_ = str; }

    Value_node( Generic_node & prev_node ):Generic_node(type_null, prev_node)                 {}
    Value_node( Generic_node & prev_node, bool b ):Generic_node(type_bool, prev_node)         { value_bool_ = b;}
    Value_node( Generic_node & prev_node, int i ):Generic_node(type_number, prev_node)        { value_number_ = i; }
    Value_node( Generic_node & prev_node, double lf ):Generic_node(type_number, prev_node)    { value_number_ = lf; }
    Value_node( Generic_node & prev_node, char * str ):Generic_node(type_string, prev_node)   { value_string_ = str; }
    Value_node( Generic_node & prev_node, string & str ):Generic_node(type_string, prev_node) { value_string_ = str; }

    virtual Generic_node & get_item( const int id ) {
        //Vector-like access to Value - return empty
        return *empty_node_generic_;
    }
    virtual Generic_node & get_key( const string & key ) {
        //Record-like access to Value - return empty
        return *empty_node_generic_;
    }
    virtual Generic_node & get_key_check( const string & key, int & err_code ) {
        //Record-like access to Value - return empty & error
        err_code = 1;
        return *empty_node_generic_;
    }
    virtual Generic_node & get_item( const size_t id, Generic_node & default_tree ) {
        //Vector-like access to Value - return empty
        return *empty_node_generic_;
    }
    virtual Generic_node & get_key( const string & key, Generic_node & default_tree ) {
        //Record-like access to Value - return empty
        return *empty_node_generic_;
    }
    virtual Generic_node & get_item_check( const size_t id, int & err_code ) {
        //Vector-like access to Value - return empty & error
        err_code = 1;
        return *empty_node_generic_;
    }

    virtual Value_node & as_value( void ) { return (*this); }
    friend ostream & operator<<( ostream & stream, Value_node & node );

    int      set_value( int i )        { value_type_ = type_number; return value_number_ = i;}
    double   set_value( double lf )    { value_type_ = type_number; return value_number_ = lf; }
    bool     set_value( bool b )       { value_type_ = type_bool;   return value_bool_   = b; }
    const char *   set_value( char * str )   { value_type_ = type_string; value_string_ = str; return value_string_.c_str(); }
    string & set_value( string & str ) { value_type_ = type_string; value_string_ = str; return value_string_; }
    void     set_null()                { value_type_ = type_null; }

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

    virtual ~Value_node();
};

}//end namespace

#endif /* VALUE_NODE_HPP_ */
