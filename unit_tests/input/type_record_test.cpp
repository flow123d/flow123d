/*
 * type_record_test.cpp
 *
 *  Created on: May 4, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>

#include <input/type_record.hh>



/**
 * Test Record class.
 */
TEST(InputTypeRecord, declare_key_scalars) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";


   // make auxiliary record and test declare_key for
   // - various Scalar types (excluding Selection): Integer, Bool, Double, String, FilePath
   // - various decalre_key templates: with/without default, shared_ptr/ reference
   // - various default values
   enum Colors {
	  white, black, red
   };

   Integer digits_type(0, 8);
   Double time_type(0.0);
   Selection sel = Selection("Color selection")
      .add_value(black, "black")
      .add_value(red, "red")
      .close();

   static Record rec = Record("SomeRecord", "desc.")
   	  .declare_key("file", FileName::output(), Default::optional(), "desc.")
      .declare_key("digits",digits_type, Default("8"), "desc.")
      .declare_key("compression", Bool(),"desc.")
      .declare_key("start_time", time_type,"desc.")
      .declare_key("data_description", String(), Default::optional(),"")
	  .declare_key("plot_color", sel, "Color to plot the fields in file.")
	  .close();


   // errors during declaration
   Record rec_empty;

#ifdef FLOW123D_DEBUG_ASSERTS
   EXPECT_THROW_WHAT( {rec_empty.declare_key("xx", Integer(), "");}, ExcXprintfMsg, ".*into closed record 'EmptyRecord'.");
#endif

   Record rec_fin = Record("xx","")
	  .close();
   EXPECT_THROW_WHAT( {rec_fin.declare_key("xx", String(),"");}, ExcXprintfMsg, "Can not add .* into closed record");


//   This no more fails: Declaration of incomplete (unfinished) keys is possible.
   /*Record rec_unfin("yy","");
   rec.declare_key("yy", rec_unfin, "");*/

   EXPECT_THROW_WHAT( { Record rec_redeclare = Record("yy","")
							.declare_key("data_description", String(), Default::optional(),"")
							.declare_key("data_description", String(),"")
							.close();
   	   	   	   	   	  }, ExcXprintfMsg, "Re-declaration of the key:");

   EXPECT_THROW_WHAT( { Record("yy","")
   	   	   	   				.declare_key("wrong_double", Double(), Default("1.23 4"),"")
							.close();
   	   	   	   	   	  }, ExcWrongDefault, "Default value .* do not match type: 'Double';");

/*
   // test documentation of default_at_read_time
   {
       Record rec("Rec", "");
       rec.declare_key("int_key", Integer(), Default::read_time("Default value provided at read time."), "");
       rec.close();

       stringstream out;
       out << rec;
       EXPECT_EQ("\nRecord 'Rec' (1 keys).\n# \n----------\n    int_key = \"Default value provided at read time.\" is Integer in [-2147483648, 2147483647]\n---------- Rec\n\n",
                 out.str());
   }

   {
       // test documentation of OPTIONAL keys
       Record rec("Rec", "");
       rec.declare_key("int_key", Integer(), Default::optional(), "Doc");
       rec.close();

       stringstream out;
       out << rec;
       EXPECT_EQ("\nRecord 'Rec' (1 keys).\n# \n----------\n    int_key = <OPTIONAL> is Integer in [-2147483648, 2147483647]\n              # Doc\n---------- Rec\n\n",
                 out.str());
   }

   {
       // test documentation of OPTIONAL keys
       Record rec("Rec", "");
       rec.declare_key("int_key", Integer(), Default::obligatory(), "Doc");
       rec.close();

       stringstream out;
       out << rec;
       EXPECT_EQ("\nRecord 'Rec' (1 keys).\n# \n----------\n    int_key = <OBLIGATORY>\n *is Integer in [-2147483648, 2147483647]\n              # Doc\n---------- Rec\n\n",
                 out.str());
   }
*/
}

TEST(InputTypeRecord, declare_key_arrays) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    // array type passed through shared_ptr
    Array array_of_int(Integer(0), 5, 100 );

    static Record array_record = Record("RecordOfArrays", "desc.")
    	.declare_key("array_of_5_ints", array_of_int,"Some bizare array.")
    	// array type passed by reference
    	.declare_key("array_of_str", Array( String() ),"Desc. of array")
    	.declare_key("array_of_str_1", Array( String() ), "Desc. of array")
    	// allow default values for an array
    	.declare_key("array_with_default", Array( Double() ), Default("3.2"), "");
    EXPECT_THROW_WHAT( { array_record.declare_key("some_key", Array( Integer() ), Default("ahoj"), ""); }, ExcWrongDefault,
                  "Default value 'ahoj' do not match type: 'Integer';"
                 );
    array_record.close();

    static Record array_record2 = Record("RecordOfArrays2", "desc.");
    EXPECT_THROW_WHAT( { array_record2.declare_key("some_key", Array( Double(), 2 ), Default("3.2"), ""); }, ExcWrongDefault,
                  "Default value '3.2' do not match type: 'array_of_Double';"
                 );
    array_record2.close();
}

TEST(InputTypeRecord, allow_convertible) {
using namespace Input::Type;
//::testing::FLAGS_gtest_death_test_style = "threadsafe";

    {
    Record sub_rec = Record( "SubRecord", "")
    				 .declare_key("default_bool", Bool(), Default("false"), "")
    				 .declare_key("optional_bool", Bool(), Default::optional(), "")
    				 .declare_key("read_time_bool", Bool(), Default::read_time(""), "")
     	 	 	 	 .declare_key("int_key", Integer(),  "")
     	 	 	 	 .allow_auto_conversion("int_key")
					 .close();
    sub_rec.finish();

    EXPECT_EQ(3, sub_rec.auto_conversion_key_iter()->key_index );
    }

    {
    static Record sub_rec = Record( "SubRecord", "")
    	.declare_key("obligatory_int", Integer(), Default::obligatory(), "")
    	.declare_key("int_key", Integer(),  "")
    	.allow_auto_conversion("int_key")
		.close();
    EXPECT_THROW_WHAT( {sub_rec.finish();}, ExcXprintfMsg,
    		"Finishing Record auto convertible from the key 'int_key', but other obligatory key: 'obligatory_int' has no default value.");
    }

    {
    static Record sub_rec = Record( "SubRecord", "")
    		.close();
    EXPECT_TRUE( sub_rec.finish() );
    EXPECT_EQ( sub_rec.end(), sub_rec.auto_conversion_key_iter() );
    }

}


TEST(InputTypeRecord, declare_key_record) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    // Test that Record has to be passed as shared_ptr
    //ASSERT_DEATH( {record_record->declare_key("sub_rec_1", Record("subrec_type", "desc") , "desc"); },
    //              "Complex type .* shared_ptr."
    //              );

    Record other_record = Record("OtherRecord","desc")
    	.close();

    static Record sub_rec = Record( "SubRecord", "")
    	.declare_key("bool_key", Bool(), Default("false"), "")
    	.declare_key("int_key", Integer(),  "")
    	.allow_auto_conversion("int_key")
    	.close();

    Record record_record = Record("RecordOfRecords", "")
    	.declare_key("sub_rec_1", other_record, "key desc");
	EXPECT_THROW_WHAT( { record_record.declare_key("sub_rec_2", other_record, Default("2.3"), "key desc"); }, ExcWrongDefault,
            "Default value '2.3' do not match type: 'OtherRecord';" );
	record_record.close();

    Record record_record2 = Record("RecordOfRecords2", "")
    	.declare_key("sub_rec_dflt", sub_rec, Default("123"), "");
    EXPECT_THROW_WHAT( { record_record2.declare_key("sub_rec_dflt2", sub_rec, Default("2.3"), ""); } , ExcWrongDefault,
            "Default value '2.3' do not match type: 'Integer';" );
    record_record2.close();

    // recursion  -  forbidden
    //record_record->declare_key("sub_rec_2", record_record, "desc.");

}

TEST(InputTypeRecord, iterating) {
    using namespace Input::Type;

    Record output_record;

    {
		Integer digits_type((int)0, (int)8);
		Double time_type(0.0);

		output_record = Record("OutputRecord", "Information about one file for field data.")
			.declare_key("file", FileName::output(), Default::optional(), "File for output stream.")
			.declare_key("digits",digits_type, Default("8"),
					"Number of digits used for output double values into text output files.")
			.declare_key("compression", Bool(), "Whether to use compression of output file.")
			.declare_key("start_time", time_type, "Simulation time of first output.")
			.declare_key("data_description", String(), "")
			.close();
		output_record.finish();
    } // delete local variables

    // methods begin() and end(), basic work with iterators
    Record::KeyIter it = output_record.begin();
    EXPECT_EQ( 0, it->key_index);
    EXPECT_EQ("file", it->key_);
    EXPECT_EQ("File for output stream.", it->description_);
    EXPECT_EQ(typeid(FileName) , typeid(*(it->type_)));
    EXPECT_EQ(FilePath::output_file, static_cast<const FileName *>( &(*it->type_) )->get_file_type() );
    it+=4;
    EXPECT_EQ( "data_description", it->key_);
    EXPECT_EQ( output_record.end(), it+1 );
    // method size()
    EXPECT_EQ(5, output_record.size());
    //method key_index
    EXPECT_EQ(0,output_record.key_index("file"));
    EXPECT_EQ(2,output_record.key_index("compression"));
    EXPECT_THROW({output_record.key_index("x_file");}, Record::ExcRecordKeyNotFound );
    // method key_iterator
    EXPECT_EQ(2,output_record.key_iterator("compression")->key_index);


}

TEST(InputTypeRecord, check_key_validity) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Record output_record("OutputRecord",
            "Information about one file for field data.");
    EXPECT_THROW_WHAT( {output_record.declare_key("a b",Bool(),"desc."); },
    		ExcXprintfMsg, "Invalid key identifier"
            );
    EXPECT_THROW_WHAT( {output_record.declare_key("AB",Bool(),"desc."); },
    		ExcXprintfMsg, "Invalid key identifier"
            );
    EXPECT_THROW_WHAT( {output_record.declare_key("%$a",Bool(),"desc."); },
    		ExcXprintfMsg, "Invalid key identifier"
            );

}

TEST(InputTypeRecord, RecordCopy) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Record output_record("OutputRecord", "");
    output_record.declare_key("file", FileName::output(), "");


    Integer digits_type((int)0, (int)8);

    Record copy_rec = output_record
    	.declare_key("digits",digits_type, "")
		.close();
    copy_rec.finish();

    EXPECT_EQ( true, output_record.is_finished());
    EXPECT_EQ( true, output_record.has_key("digits") );
}


TEST(InputTypeRecord, copy_keys) {
using namespace Input::Type;

    Record rec1 =
    		Record("Rec1", "")
    		.declare_key("a", Integer(), "a from rec1")
			.close();

    Record rec2 =
    		Record("Rec2","")
    		.declare_key("a", Integer(), "a from rec2")
    		.declare_key("b", Integer(), "b from rec2")
    		.declare_key("c", Integer(), "c from rec2")
			.close();

    Record composite =
    		Record("composite","")
    		.declare_key("b", Integer(), "b from composite")
    		.copy_keys(rec1)
    		.copy_keys(rec2)
			.close();

    composite.finish();

    EXPECT_FALSE(rec1.is_finished());
    EXPECT_FALSE(rec2.is_finished());

    EXPECT_EQ(3, composite.size());
    EXPECT_EQ("a from rec1", composite.key_iterator("a")->description_);
    EXPECT_EQ("b from composite", composite.key_iterator("b")->description_);
    EXPECT_EQ("c from rec2", composite.key_iterator("c")->description_);
}

