#ifndef KEYCONNECTOR_HPP_
#define KEYCONNECTOR_HPP_

#include <string>
#include "Generic_node.hpp"

using namespace std;

namespace flow {

/*!
 * @brief Key connector - used for declaration, key description and default value.
 *          It is inserted between JSON key and value in record type.
 *          Until user calls DECLARE_* function, it is not possible to read stored value.
 *
 */
class Key_connector {
    bool declared_ ;
    Generic_node * default_node_;
    string description_;
public:
    Key_connector()
    {
        declared_ = false;
        default_node_ = NULL;
    }
    Key_connector( string & desc ) {
        declared_ = true;
        description_ = desc;
        default_node_ = NULL;
    }
    Key_connector( string & desc, Generic_node & default_node ) {
        declared_ = true;
        description_ = desc;
        default_node_ = &default_node;
    }
    //TODO: Deep copy constructor. Bude vubec potreba? Je tu pointer, takze bacha...
    //zatim delam shallow kopii...

    bool is_declared() {
        return declared_;
    }

    string & get_description() {
        return description_;
    }

    Generic_node & get_default() {
        return * default_node_;
    }

    void set_declared(bool declared) {
        declared_ = declared;
    }

    void set_description(string & description) {
        description_ = description;
    }
    virtual ~Key_connector();
};


} /* namespace flow */
#endif /* KEYCONNECTOR_HPP_ */
