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


#include <gtest_throw_what.hh>
#include <vector>

#include <input/input_interface.hh>
#include <input/input_type.hh>
#include <input/type_record.hh>

#include "system/file_path.hh"

class InputInterfaceTest : public testing::Test {


protected:

    virtual void SetUp() {
        using namespace Input::Type;

        FilePath::set_io_dirs("/json_root_dir","variant_input","/output_root");

        abstr_rec_ptr = new  AbstractRecord("AbstractRecord", "desc");
        abstr_rec_ptr->finish();

        desc_a_ptr = new Record("DescendantA","");
        desc_a_ptr->derive_from(*abstr_rec_ptr);
        desc_a_ptr->declare_key("some_int", Integer(),"");
        desc_a_ptr->finish();

        desc_b_ptr = new Record("DescendantB","");
        desc_b_ptr->derive_from(*abstr_rec_ptr);
        desc_b_ptr->declare_key("some_int", Integer(),"");
        desc_b_ptr->declare_key("some_double", Double(),"");
        desc_b_ptr->finish();

        abstr_rec_ptr->no_more_descendants();

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
        main->declare_key("abstr_rec_1",*abstr_rec_ptr, "desc");
        main->declare_key("abstr_rec_2",*abstr_rec_ptr, "desc");
        main->declare_key("file_output", FileName::output(), "description");
        main->declare_key("file_input", FileName::input(), "description");
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

            StorageArray * desc_a = new StorageArray(2);
            desc_a->new_item(0,new StorageInt(0));
            desc_a->new_item(1,new StorageInt(234));

            StorageArray * desc_b = new StorageArray(3);
            desc_b->new_item(0,new StorageInt(1));
            desc_b->new_item(1,new StorageInt(345));
            desc_b->new_item(2,new StorageDouble(3.45));


            StorageArray * main_array = new StorageArray(11);
            main_array->new_item(0, sub_rec->deep_copy());
            main_array->new_item(1, sub_array_int->deep_copy());
            main_array->new_item(2, sub_array_sub_rec->deep_copy());
            main_array->new_item(3, new StorageInt(456));
            main_array->new_item(4, new StorageDouble(4.56));
            main_array->new_item(5, new StorageBool(false));
            main_array->new_item(6, new StorageString("456"));
            main_array->new_item(7, desc_a);
            main_array->new_item(8, desc_b);
            main_array->new_item(9, new StorageString("output_subdir/output.vtk"));
            main_array->new_item(10, new StorageString("input/${INPUT}/input_subdir/input.in"));

            delete sub_array_int;
            delete sub_rec;
            delete sub_array_sub_rec;

            this->storage = main_array;

        }

    }

    virtual void TearDown() {
        delete main;
        delete storage;
        delete desc_a_ptr;
        delete desc_b_ptr;
        delete abstr_rec_ptr;
    };

    ::Input::Type::Record *main;
    ::Input::Storage * storage;

    ::Input::Type::Record *desc_a_ptr;
    ::Input::Type::Record *desc_b_ptr;
    ::Input::Type::AbstractRecord *abstr_rec_ptr;
};

TEST_F(InputInterfaceTest, ReadFromRecord) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    using namespace Input;

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

    EXPECT_FALSE( record.key<bool>("some_bool") );

    EXPECT_EQ("456", record.key<string>("some_string") );

    EXPECT_EQ("/output_root/output_subdir/output.vtk", (string) record.key<FilePath>("file_output") );
    EXPECT_EQ("/json_root_dir/input/variant_input/input_subdir/input.in", (string) record.key<FilePath>("file_input") );


    // read record
    Record sub_record( record.key<Record>("some_record") );
        i = sub_record.key<int>("some_integer");
    EXPECT_EQ(123, i);

    // read array key
    Array array = record.key<Array>("array_of_int");
    EXPECT_EQ(2, *( ++array.begin<int>()) );

    EXPECT_THROW_WHAT( { record.key<int>("unknown_key");} , Input::Type::Record::ExcRecordKeyNotFound, "Key 'unknown_key' not found in Record");
    EXPECT_THROW_WHAT( { record.key<int>("some_bool");} , Input::ExcTypeMismatch , "can not make iterator with type");


}

struct Data {
    bool b;
    int i;
    double d;
    string s;
};

TEST_F(InputInterfaceTest, ReadFromArray) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    using namespace Input;

    Record record(storage, *main);
    Array array = record.key<Array>("array_of_int");

    std::vector<int> vec_int;
    // reading scalars form array - manually
    for(Iterator<int> it=array.begin<int>(); it != array.end(); ++it) vec_int.push_back(*it);
    EXPECT_EQ(2, vec_int.size());
    EXPECT_EQ(1, vec_int[0]);
    EXPECT_EQ(2, vec_int[1]);

    // using copy_to method template
    array.copy_to(vec_int);
    EXPECT_EQ(2, vec_int.size());
    EXPECT_EQ(1, vec_int[0]);
    EXPECT_EQ(2, vec_int[1]);

    Iterator<int> it = array.begin<int>();
    ++it;
    ++it;
    ++it;
    EXPECT_DEATH( {int ii = *it;}, "out of array of size:");




    array = record.key<Array>("array_of_sub_rec");
    Data * data_array = new Data[array.size()];

    // can not use copy_to for this type !!!
    int idx=0;
    for(Iterator<Record> it = array.begin<Record>(); it != array.end(); ++it, ++idx ) {
        data_array[idx].b = it->key<bool>("some_bool");
        if (it->has_key("some_int", data_array[idx].i) ) {
            EXPECT_EQ(123,data_array[idx].i);
        }
        it->has_key("some_double", data_array[idx].d);
        EXPECT_EQ(1.23, data_array[idx].d);
        it->has_key("some_string", data_array[idx].s);
        EXPECT_EQ("123", data_array[idx].s);
    }


    delete[] data_array;
}

TEST_F(InputInterfaceTest, ReadFromAbstract) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    using namespace Input;

    Record record(storage, *main);

    {
        AbstractRecord a_rec = record.key<AbstractRecord>("abstr_rec_1");
        EXPECT_EQ(a_rec.type(), *( this->desc_a_ptr ));
        if (a_rec.type() == *( this->desc_a_ptr )) {
            Record rec( a_rec );
            EXPECT_EQ(234, rec.key<int>("some_int") );
        } else
        if (a_rec.type() == *( this->desc_b_ptr )) {
            Record rec( a_rec );
            EXPECT_EQ(345, rec.key<int>("some_int") );
            EXPECT_EQ(3.45, rec.key<int>("some_double") );
        }

    }


    {
        AbstractRecord a_rec = record.key<AbstractRecord>("abstr_rec_2");
        EXPECT_EQ(a_rec.type(), *( this->desc_b_ptr ));
        if (a_rec.type() == *( this->desc_a_ptr )) {
            Record rec( a_rec );
            EXPECT_EQ(234, rec.key<int>("some_int") );
        } else
        if (a_rec.type() == *( this->desc_b_ptr )) {
            Record rec( a_rec );
            EXPECT_EQ(345, rec.key<int>("some_int") );
            EXPECT_EQ(3.45, rec.key<double>("some_double") );
        }

    }
}


