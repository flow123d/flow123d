/*
 * input_interface_test.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *
 * TODO:
 *
 *  - klice nezadavat jako retezce ale pres makro KEY(nazev_klice) to
 *    umozni pocitani hashu pri kompilaci pripadne dalsi kejkle
 *  - implementovat deklaraci Selection (pouziti enum_macro zatim nevhodne, protoze to je nestabilni, leda casem vytvorit neco vlastniho)
 *    zatim rucne vkladat enum klice do Selection/ AbstracRecord
 *  - implementovat cteni key<enum ..>
 *  - implementovat deklaraci AbstractRecord a jeho cteni jak key<Record> tak key<enum...>
 *  - implementovat Iterator<Record/Array>::operator->
 *
 *
 *  - nevytvaret deklaraci Recordu a dalsich typu runtime, ale pri kompilaci
 *    tedy jako skutecnou hierarchii trid, to by umoznilo statickou kontrolu
 *    kompatibility typu v konstruktoru Iterator<T>
 *
 */


#include <gtest/gtest.h>
#include <vector>

#include <input/input_interface.hh>
#include <input/input_type.hh>

class InputInterfaceTest : public testing::Test {


protected:

    virtual void SetUp() {
        using namespace Input::Type;

        // declare structure of input file
        main = boost::make_shared<Record>("MainRecord","desc");

        boost::shared_ptr<Record> sub_record = boost::make_shared<Record>("SubRecord","desc");
        sub_record->declare_key("array_of_int", Array(Integer()), "desc");
        sub_record->declare_key("some_integer", Integer(), "desc");
        sub_record->declare_key("some_double", Double(), "desc");
        sub_record->declare_key("some_bool", Bool(), "desc");
        sub_record->declare_key("some_string", String(), "desc");

        main->declare_key("some_record", sub_record, "desc");
        main->declare_key("array_of_int", Array(Integer()), "desc");
        main->declare_key("array_of_sub_rec", Array( sub_record ), "desc");
        main->declare_key("some_integer", Integer(), "desc");
        main->declare_key("some_double", Double(), "desc");
        main->declare_key("some_bool", Bool(), "desc");
        main->declare_key("some_string", String(), "desc");

        // read data storage tree
        this->storage = new Input::Interface::Storage();

    }

    virtual void TearDown() {
        main.reset(); // effectively lose all declarations
    };

    boost::shared_ptr< ::Input::Type::Record > main;
    ::Input::Interface::Storage * storage;
};

TEST_F(InputInterfaceTest, test_reading) {
    using namespace Input::Interface;

    Record record(*storage, *main);

    // read scalar keys
    int i;
        i = record.key<int>("some_integer");
        i = record.key<char>("some_integer");
        i = record.key<short int>("some_integer");

    double d = record.key<double>("some_double");
           d = record.key<float>("some_double");

    bool b = record.key<bool>("some_bool");
    string s = record.key<string>("some_string");

    // read record
    Record sub_record( record.key<Record>("some_record") );
        i = sub_record.key<int>("some_integer");

    // read array key
    Array array = record.key<Array>("array_of_int");

    // reading scalars form array
    std::vector<int> vec_int;
    //array.copy_to(vec_int);

    for(Iterator<int> it=array.begin<int>(); it != array.end(); ++it) vec_int.push_back(*it);

    array = record.key<Array>("array_of_sub_rec");

    // can not use copy_to for this type !!!
    for(Iterator<Record> it = array.begin<Record>(); it != array.end(); ++ it ) {
        i = (*it).key<int>("some_integer");
    }



}
