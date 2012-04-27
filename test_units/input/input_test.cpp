#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cfloat>

#include "../input.hpp" //The only necessary header file.

using namespace std;
using namespace flow;

int main()
{
    //non-existent node example

    Generic_node gnode;
    cout << endl << "Default demo:" << endl;
    int default_int = 12321;
    int returned_int;
    returned_int = gnode.get_key("foo").get_item(10).as_value().get_int(default_int);
    cout << "default_int=" << default_int << endl;
    cout << "returned_int=" << returned_int << endl;

    //access to instance as an ancestor

    Value_node vnode;
    Generic_node & gnode_r = vnode;
    Generic_node * gnode_p;
    gnode_p = new Value_node(123);
    cout << endl << "Access as ancestor demo:" << endl;
    cout << "gnode " << gnode << endl;
    cout << "vnode " << vnode << endl;
    cout << "gnode reference to vnode " << gnode_r << endl;
    cout << "gnode pointer to vnode " << (*gnode_p) << endl;

    return 0;
}
