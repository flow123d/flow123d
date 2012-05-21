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

        {
        // declare structure of input file
        this->main = new Record("MainRecord","desc");

        Record sub_record("SubRecord","desc");
        sub_record.declare_key("array_of_int", Array(Integer()), "desc");
        sub_record.declare_key("some_integer", Integer(), "desc");
        sub_record.declare_key("some_double", Double(), "desc");
        sub_record.declare_key("some_bool", Bool(), "desc");
        sub_record.declare_key("some_string", String(), "desc");
        sub_record.finish();

        main->declare_key("some_record", sub_record, "desc");
        main->declare_key("array_of_int", Array(Integer()), "desc");
        main->declare_key("array_of_sub_rec", Array( sub_record ), "desc");
        main->declare_key("some_integer", Integer(), "desc");
        main->declare_key("some_double", Double(), "desc");
        main->declare_key("some_bool", Bool(), "desc");
        main->declare_key("some_string", String(), "desc");
        main->finish();
        }

        // construct some storage


        {
            using namespace Input;

            StorageArray * sub_array_int = new StorageArray(2);
            sub_array_int->new_item(0, new StorageInt(1));
            sub_array_int->new_item(1, new StorageInt(2));

            StorageArray * sub_rec = new StorageArray(5);
            sub_rec->new_item(0, sub_array_int->deep_copy());
            sub_rec->new_item(1, new StorageInt(123));
            sub_rec->new_item(2, new StorageDouble(1.23));
            sub_rec->new_item(3, new StorageBool(true));
            sub_rec->new_item(4, new StorageString("123"));

            StorageArray * sub_array_sub_rec = new StorageArray(2);
            sub_array_sub_rec->new_item(0, sub_rec->deep_copy());
            sub_array_sub_rec->new_item(1, sub_rec->deep_copy());

            StorageArray * main_array = new StorageArray(7);
            main_array->new_item(0, sub_rec->deep_copy());
            main_array->new_item(1, sub_array_int->deep_copy());
            main_array->new_item(2, sub_array_sub_rec->deep_copy());
            main_array->new_item(3, new StorageInt(456));
            main_array->new_item(4, new StorageDouble(4.56));
            main_array->new_item(5, new StorageBool(false));
            main_array->new_item(6, new StorageString("456"));

            delete sub_array_int;
            delete sub_rec;
            delete sub_array_sub_rec;

            this->storage = main_array;

        }

    }

    virtual void TearDown() {
        delete main;
        delete storage;
    };

    ::Input::Type::Record *main;
    ::Input::Interface::Storage * storage;
};

TEST_F(InputInterfaceTest, ReadFromRecord) {
    using namespace Input::Interface;

    Record record(storage, *main);

    // read scalar keys
    int i;
        i = record.key<int>("some_integer");
        EXPECT_EQ(456,i);
        i = record.key<char>("some_integer");
        EXPECT_EQ((char)(456),i);
        i = record.key<short int>("some_integer");
        EXPECT_EQ((short int)(456),i);

    double d = record.key<double>("some_double");
           EXPECT_EQ(4.56, d);
           d = record.key<float>("some_double");
           EXPECT_EQ((float)4.56, d);

    EXPECT_EQ(true, record.has_key("some_double", d));
    EXPECT_EQ(4.56, d);

    bool b = record.key<bool>("some_bool");
    string s = record.key<string>("some_string");

    // read record
    Record sub_record( record.key<Record>("some_record") );
        i = sub_record.key<int>("some_integer");

    // read array key
    Array array = record.key<Array>("array_of_int");

}

struct Data {
    bool b;
    int i;
    double d;
    string s;
};

TEST_F(InputInterfaceTest, ReadFromArray) {
    using namespace Input::Interface;

    Record record(storage, *main);
    Array array = record.key<Array>("array_of_int");

    std::vector<int> vec_int;
    // reading scalars form array - manually
    for(Iterator<int> it=array.begin<int>(); it != array.end(); ++it) vec_int.push_back(*it);
    // using copy_to method template
    array.copy_to(vec_int);



    array = record.key<Array>("array_of_sub_rec");

    Data * data_array = new Data[array.size()];

    // can not use copy_to for this type !!!
    int idx=0;
    for(Iterator<Record> it = array.begin<Record>(); it != array.end(); ++it, ++idx ) {
        data_array[idx].b = it->key<bool>("some_bool");
        if (it->has_key("some_int", data_array[idx].i) ) {
            xprintf(Msg, "OK.");
        }
        it->has_key("some_double", data_array[idx].d);
        it->has_key("some_string", data_array[idx].s);

    }

    delete[] data_array;
}

