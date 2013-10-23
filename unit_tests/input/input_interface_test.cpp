/*
 * input_interface_test.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */


#include <gtest_throw_what.hh>
#include <vector>

#include <input/accessors.hh>
#include <input/input_type.hh>
#include <input/type_record.hh>

#include "system/file_path.hh"

    enum SelectionToRead {
        value_a = 0,
        value_b = 1,
        value_c = 2
    };


class InputInterfaceTest : public testing::Test {
public:


protected:

    virtual void SetUp() {
        using namespace Input::Type;

        FilePath::set_io_dirs("/json_root_dir","/json_root_dir","variant_input","/output_root");

        abstr_rec_ptr = new  AbstractRecord("AbstractRecord", "desc");
        abstr_rec_ptr->finish();

        selection_ptr = new Selection("NameOfSelectionType");
        selection_ptr->add_value(value_a, "A", "");
        selection_ptr->add_value(value_b, "B", "");
        selection_ptr->add_value(value_c, "C", "");
        selection_ptr->finish();

        desc_a_ptr = new Record("DescendantA","");
        desc_a_ptr->derive_from(*abstr_rec_ptr);
        desc_a_ptr->declare_key("some_int", Integer(),Default("1"),"");
        desc_a_ptr->finish();

        desc_b_ptr = new Record("DescendantB","");
        desc_b_ptr->derive_from(*abstr_rec_ptr);
        desc_b_ptr->declare_key("some_int", Integer(),Default("2"),"");
        desc_b_ptr->declare_key("some_double", Double(),Default::obligatory(),"");
        desc_b_ptr->finish();

        abstr_rec_ptr->no_more_descendants();

        {
        // declare structure of input file
        this->main = new Record("MainRecord","desc");

        Record sub_record("SubRecord","desc");
        sub_record.declare_key("array_of_int", Array(Integer()), "desc");
        sub_record.declare_key("some_integer", Integer(), Default::obligatory(), "desc");
        sub_record.declare_key("some_double", Double(), "desc");
        sub_record.declare_key("some_bool", Bool(), Default("true"), "desc");
        sub_record.declare_key("some_string", String(), "desc");
        sub_record.finish();



        main->declare_key("some_record", sub_record, Default::obligatory(), "desc");
        main->declare_key("array_of_int", Array(Integer()),Default::obligatory(), "desc");
        main->declare_key("array_of_sub_rec", Array( sub_record ),Default::obligatory(), "desc");
        main->declare_key("some_integer", Integer(),Default::obligatory(), "desc");
        main->declare_key("some_double", Double(),Default::obligatory(), "desc");
        main->declare_key("some_bool", Bool(),Default::obligatory(), "desc");
        main->declare_key("some_string", String(),Default::obligatory(), "desc");
        main->declare_key("abstr_rec_1",*abstr_rec_ptr,Default::obligatory(), "desc");
        main->declare_key("abstr_rec_2",*abstr_rec_ptr,Default::obligatory(), "desc");
        main->declare_key("file_output", FileName::output(),Default::obligatory(), "description");
        main->declare_key("file_input", FileName::input(),Default::obligatory(), "description");
        main->declare_key("optional_int", Integer(), "");
        main->declare_key("selection", *selection_ptr, Default::obligatory(), "");
        main->declare_key("default_int", Integer(), Default("1234"), "");
        main->declare_key("optional_int2", Integer(), "");
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


            StorageArray * main_array = new StorageArray(15);
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
            main_array->new_item(11, new StorageNull());
            main_array->new_item(12, new StorageInt(1));
            main_array->new_item(13, new StorageInt(1234));
            main_array->new_item(14, new StorageInt(12345));

            // check copy constructors and pimpl implementation of Record
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
        delete selection_ptr;
    };

    ::Input::Type::Record *main;
    ::Input::StorageBase * storage;

    ::Input::Type::Record *desc_a_ptr;
    ::Input::Type::Record *desc_b_ptr;
    ::Input::Type::AbstractRecord *abstr_rec_ptr;

    ::Input::Type::Selection *selection_ptr;
};

TEST_F(InputInterfaceTest, RecordVal) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    using namespace Input;

    Address addr(storage, main);
    Record record(addr, *main);

    // read scalar keys
    int i;
        i = record.val<int>("some_integer");
        EXPECT_EQ(456,i);
        i = record.val<char>("some_integer");
        EXPECT_EQ((char)(456),i);
        i = record.val<short int>("some_integer");
        EXPECT_EQ((short int)(456),i);

    double d = record.val<double>("some_double");
           EXPECT_EQ(4.56, d);
           d = record.val<float>("some_double");
           EXPECT_EQ((float)4.56, d);

    //EXPECT_EQ(true, record.has_key("some_double", d));
    //EXPECT_EQ(4.56, d);

    EXPECT_FALSE( record.val<bool>("some_bool") );

    EXPECT_EQ("456", record.val<string>("some_string") );

    EXPECT_EQ("/output_root/output_subdir/output.vtk", (string) record.val<FilePath>("file_output") );
    EXPECT_EQ("/json_root_dir/input/variant_input/input_subdir/input.in", (string) record.val<FilePath>("file_input") );

    // read enum from selection
    EXPECT_EQ( value_b, record.val<SelectionToRead>("selection") );

    // read record
    Record sub_record( record.val<Record>("some_record") );
        i = sub_record.val<int>("some_integer");
    EXPECT_EQ(123, i);

    // read array key
    Array array = record.val<Array>("array_of_int");
    EXPECT_EQ(2, *( ++array.begin<int>()) );

    EXPECT_THROW_WHAT( {record.val<string>("some_double");}, ExcTypeMismatch,
            "Program Error: Key:'some_double'. Can not construct Iterator<T> with C.. type T='Ss';");
    EXPECT_THROW( {record.val<string>("unknown");}, Type::Record::ExcRecordKeyNotFound );

#ifdef DEBUG_ASSERTS
    EXPECT_THROW_WHAT( {record.val<int>("optional_int");}, ExcAssertMsg,
            "The key 'optional_int' is declared as optional .*you have to use Record::find instead.");
#endif

}

TEST_F(InputInterfaceTest, RecordFind) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    using namespace Input;

    Address addr(storage, main);
    Record record(addr, *main);

    Address addr_child_rec( *(addr.down(0)) );
    const Input::Type::Record * child_type_rec = static_cast<const Type::Record *>( main->begin()->type_.get() );
    Record record_child_rec(addr_child_rec, *child_type_rec );

    Address addr_child_int( *(record_child_rec.get_address().down(1)) );

    EXPECT_EQ("/some_record", record_child_rec.get_address().make_full_address());
    EXPECT_EQ("/some_record/some_integer", addr_child_int.make_full_address());

    // read scalar keys

    Iterator<int> it = record.find<int>("some_integer");
    EXPECT_TRUE(it);
    EXPECT_EQ(456, *it);

    it = record.find<int>("optional_int");
    EXPECT_FALSE(it);

    Iterator<Record> it_r = record.find<Record>("some_record");
    EXPECT_EQ(123, it_r->val<int>("some_integer"));

    EXPECT_THROW_WHAT( {record.find<string>("some_double");}, ExcTypeMismatch,
            "Program Error: Key:'some_double'. Can not construct Iterator<T> with C.. type T='Ss';");
    EXPECT_THROW( {record.find<string>("unknown");}, Type::Record::ExcRecordKeyNotFound );

}

TEST_F(InputInterfaceTest, Record_opt_val) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    using namespace Input;

    Address addr(storage, main);
    Record record(addr, *main);

    // read obligatory or default key
    int value;
    EXPECT_TRUE(record.opt_val("some_integer", value));
    EXPECT_EQ(456, value);
    EXPECT_TRUE(record.opt_val("default_int", value));
    EXPECT_EQ(1234, value);

    // read optional key
    EXPECT_FALSE(record.opt_val("optional_int", value));
    EXPECT_TRUE(record.opt_val("optional_int2", value));
    EXPECT_EQ(12345, value);
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

    Address addr(storage, main);
    Record record(addr, *main);
    Array array = record.val<Array>("array_of_int");

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




    array = record.val<Array>("array_of_sub_rec");
    Data * data_array = new Data[array.size()];

    // can not use copy_to for this type !!!
    int idx=0;
    for(Iterator<Record> it = array.begin<Record>(); it != array.end(); ++it, ++idx ) {
        data_array[idx].b = it->val<bool>("some_bool");

//        if (it->has_key("some_int", data_array[idx].i) ) {
//            EXPECT_EQ(123,data_array[idx].i);
 //       }
 //       it->has_key("some_double", data_array[idx].d);
 //       EXPECT_EQ(1.23, data_array[idx].d);
 //       it->has_key("some_string", data_array[idx].s);
 //       EXPECT_EQ("123", data_array[idx].s);
    }

    // check creation of empty accessor and defautl iterator
    Array empty_array;
    EXPECT_EQ(0, empty_array.size());

    Iterator<Record> it_r;


    delete[] data_array;
}

TEST_F(InputInterfaceTest, ReadFromAbstract) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    using namespace Input;

    Address addr(storage, main);
    Record record(addr, *main);

    {
        AbstractRecord a_rec = record.val<AbstractRecord>("abstr_rec_1");
        EXPECT_EQ(a_rec.type(), *( this->desc_a_ptr ));
        if (a_rec.type() == *( this->desc_a_ptr )) {
            Record rec( a_rec );
            EXPECT_EQ(234, rec.val<int>("some_int") );
        } else
        if (a_rec.type() == *( this->desc_b_ptr )) {
            Record rec( a_rec );
            EXPECT_EQ(345, rec.val<int>("some_int") );
            EXPECT_EQ(3.45, rec.val<int>("some_double") );
        }

    }


    {
        AbstractRecord a_rec = record.val<AbstractRecord>("abstr_rec_2");
        EXPECT_EQ(a_rec.type(), *( this->desc_b_ptr ));
        if (a_rec.type() == *( this->desc_a_ptr )) {
            Record rec( a_rec );
            EXPECT_EQ(234, rec.val<int>("some_int") );
        } else
        if (a_rec.type() == *( this->desc_b_ptr )) {
            Record rec( a_rec );
            EXPECT_EQ(345, rec.val<int>("some_int") );
            EXPECT_EQ(3.45, rec.val<double>("some_double") );
        }

    }
}


