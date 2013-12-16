/*
 * storage_test.cpp
 *
 *  Created on: May 5, 2012
 *      Author: jb
 */

#include <flow_gtest.hh>

#include "input/storage.hh"

//use namespace std;

TEST(Storage, all) {
using namespace Input;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    StorageArray array(7);
    array.new_item(0, new StorageNull());
    array.new_item(1, new StorageBool(true));
    array.new_item(2, new StorageInt(123));
    array.new_item(3, new StorageDouble(3.14));
    array.new_item(4, new StorageString("ahoj"));

    StorageArray * sub_array1= new StorageArray(2);
    array.new_item(5, sub_array1);
    sub_array1->new_item(0, new StorageNull());

    sub_array1=new StorageArray(2);
    sub_array1->new_item(0, new StorageInt(321));
    sub_array1->new_item(1, new StorageInt(231));
    array.new_item(6, sub_array1);

#ifdef DEBUG_ASSERTS
    EXPECT_DEATH( {array.new_item(7, sub_array1);}, "out of array of size:");
#endif

    EXPECT_TRUE(array.get_item(0)->is_null());
    EXPECT_FALSE(array.get_item(1)->is_null());
    EXPECT_EQ(true, array.get_item(1)->get_bool());
    EXPECT_EQ(123, array.get_item(2)->get_int());
    EXPECT_EQ(3.14, array.get_item(3)->get_double());
    EXPECT_EQ("ahoj", array.get_item(4)->get_string());

    EXPECT_EQ(7,array.get_array_size());

    EXPECT_TRUE(array.get_item(5)->get_item(0)->is_null());
    //EXPECT_EQ(NULL, array.get_item(5)->get_item(1) );
    EXPECT_EQ(321, array.get_item(6)->get_item(0)->get_int());
    EXPECT_EQ(231, array.get_item(6)->get_item(1)->get_int());

    EXPECT_DEATH( {array.get_item(7);} , "out of array of size:");

    // StorageArray
    EXPECT_THROW( {array.get_int();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_bool();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_double();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_string();}, ExcStorageTypeMismatch);

    // StorageNull
    EXPECT_THROW( {array.get_item(0)->get_int();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(0)->get_bool();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(0)->get_double();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(0)->get_string();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(0)->get_item(0);}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(0)->get_array_size();}, ExcStorageTypeMismatch);

    //StorageBool
     EXPECT_THROW( {array.get_item(1)->get_int();}, ExcStorageTypeMismatch);
     EXPECT_THROW( {array.get_item(1)->get_double();}, ExcStorageTypeMismatch);
     EXPECT_THROW( {array.get_item(1)->get_string();}, ExcStorageTypeMismatch);
     EXPECT_THROW( {array.get_item(1)->get_item(0);}, ExcStorageTypeMismatch);
     EXPECT_THROW( {array.get_item(1)->get_array_size();}, ExcStorageTypeMismatch);

    //StorageInt
    EXPECT_THROW( {array.get_item(2)->get_bool();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(2)->get_double();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(2)->get_string();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(2)->get_item(0);}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(2)->get_array_size();}, ExcStorageTypeMismatch);

    // StorageDouble
    EXPECT_THROW( {array.get_item(3)->get_int();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(3)->get_bool();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(3)->get_string();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(3)->get_item(0);}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(3)->get_array_size();}, ExcStorageTypeMismatch);

    // StorageString
    EXPECT_THROW( {array.get_item(4)->get_int();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(4)->get_bool();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(4)->get_double();}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(4)->get_item(0);}, ExcStorageTypeMismatch);
    EXPECT_THROW( {array.get_item(4)->get_array_size();}, ExcStorageTypeMismatch);
}

