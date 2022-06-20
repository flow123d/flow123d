/*
 * type_record_test.cpp
 *
 *  Created on: May 4, 2012
 *      Author: jb
 */


#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>


#include <input/type_record.hh>
#include <input/reader_internal_base.hh>

// Test of correct includes in type_record.hh
TEST(InputTypeRecord, includes) {
	using namespace Input::Type;
	Record rec = Record("EmptyRec", "description").close();
	EXPECT_EQ( rec.class_name(), "Record");
}


#include <input/input_type.hh>
#include <input/reader_to_storage.hh>



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
#ifdef FLOW123D_DEBUG_ASSERTS
   Record rec_empty;
   EXPECT_ASSERT_DEATH( {rec_empty.declare_key("xx", Integer(), "");}, "key : 'xx'");

   Record rec_fin = Record("xx","").close();
   EXPECT_ASSERT_DEATH( {rec_fin.declare_key("xx", String(),"");}, "key : 'xx'");
#endif


//   This no more fails: Declaration of incomplete (unfinished) keys is possible.
   /*Record rec_unfin("yy","");
   rec.declare_key("yy", rec_unfin, "");*/

   EXPECT_THROW_WHAT( { Record rec_redeclare = Record("yy","")
							.declare_key("data_description", String(), Default::optional(),"")
							.declare_key("data_description", String(),"")
							.close();
   	   	   	   	   	  }, feal::Exc_assert, "Re-declaration of the key");

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
    	.declare_key("array_with_default", Array( Double() ), Default("3.2"), "")
    	.close();
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

    EXPECT_EQ(4, sub_rec.auto_conversion_key_iter()->key_index );
    }

    {
    static Record sub_rec = Record( "SubRecord", "")
    	.declare_key("obligatory_int", Integer(), Default::obligatory(), "")
    	.declare_key("int_key", Integer(),  "")
    	.allow_auto_conversion("int_key")
		.close();
    EXPECT_THROW_WHAT( {sub_rec.finish();}, feal::Exc_assert,
    		"other_key : 'obligatory_int'");
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

    Record sub_rec = Record( "SubRecord", "")
    	.declare_key("bool_key", Bool(), Default("false"), "")
    	.declare_key("int_key", Integer(),  "")
    	.allow_auto_conversion("int_key")
    	.close();
    sub_rec.finish();

    Record record_record2 = Record("RecordOfRecords2", "")
    	.declare_key("sub_rec_dflt", sub_rec, Default("123"), "")
		.declare_key("sub_rec_dflt2", sub_rec, Default("2.3"), "")
    	.close();
    EXPECT_THROW_WHAT( { record_record2.finish(); } , Input::ReaderInternalBase::ExcAutomaticConversionError,
            "Error during automatic conversion of SubRecord record" );


    // recursion  -  forbidden
    //record_record->declare_key("sub_rec_2", record_record, "desc.");

}

TEST(InputTypeRecord, check_default_validity) {
using namespace Input::Type;

	{ // wrong default value of double
		Record rec = Record("yy","")
			.declare_key("wrong_double", Double(), Default("1.23 4"),"")
			.close();
		EXPECT_THROW_WHAT( { rec.finish(); }, ExcWrongDefaultJSON,
			"Not valid JSON of Default value '1.23 4' of type 'Double';");
	}

	{ // wrong default value of array
		Record rec = Record("SomeRecordOfArrays", "desc.")
			.declare_key("array_of_str", Array( String() ),"Desc. of array")
			.declare_key("some_key", Array( Integer() ), Default("ahoj"), "")
			.close();
		EXPECT_THROW_WHAT( { rec.finish(); }, ExcWrongDefaultJSON,
			"Not valid JSON of Default value 'ahoj' of type 'array_of_Integer';");
	}

	{ // wrong default value of array with minimal size 2
		Record rec = Record("RecordOfDoubleArray", "desc.")
			.declare_key("some_key", Array( Double(), 2 ), Default("3.2"), "")
			.close();
		EXPECT_THROW_WHAT( { rec.finish(); }, ExcWrongDefault,
					  "Default value '3.2' do not match type: 'array_of_Double';");
	}

	{ // wrong default value of SubRecord
		Record other_record = Record("OtherRecord","desc")
			.close();
		other_record.finish();

		Record rec = Record("RecordOfRecords", "")
			.declare_key("sub_rec_1", other_record, "key desc")
			.declare_key("sub_rec_2", other_record, Default("2.3"), "key desc")
			.close();
		EXPECT_THROW_WHAT( { rec.finish(); }, ExcWrongDefault,
				"Default value '2.3' do not match type: 'OtherRecord';" );
	}

	{ // generic record - passed default value
		Record generic_rec = Record("VariableRecord1","")
			.root_of_generic_subtree()
			.declare_key("key", Array(Parameter("param")), Default("[]"),"")
			.close();

		std::vector<TypeBase::ParameterPair> param_vec;
		param_vec.push_back( std::make_pair("param", std::make_shared<Integer>()) );
		Record rec = Record("RecordWithGeneric1", "")
			.declare_key("int_array", Instance(generic_rec, param_vec).close(), "key desc")
			.close();

		// simulate order and is_finished parameter of calling finish methods of IST
		generic_rec.finish(FinishStatus::generic_);
		rec.finish();
	}

	{ // generic record - wrong size of default array
		Record generic_rec = Record("VariableRecord2","")
			.root_of_generic_subtree()
			.declare_key("key", Array(Parameter("param"), 2), Default("[]"),"")
			.close();

		std::vector<TypeBase::ParameterPair> param_vec;
		param_vec.push_back( std::make_pair("param", std::make_shared<Integer>()) );
		Record rec = Record("RecordWithGeneric2", "")
			.declare_key("int_array2", Instance(generic_rec, param_vec).close(), "key desc")
			.close();

		generic_rec.finish(FinishStatus::generic_);
		EXPECT_THROW_WHAT( { rec.finish(); }, ExcWrongDefault,
				"do not match type: 'array_of_Integer';" );
	}

	{ // generic record - wrong type of default array
		Record generic_rec = Record("VariableRecord3","")
			.root_of_generic_subtree()
			.declare_key("key", Array(Parameter("param")), Default("[ 0.5 ]"),"")
			.close();

		std::vector<TypeBase::ParameterPair> param_vec;
		param_vec.push_back( std::make_pair("param", std::make_shared<Integer>()) );
		Record rec = Record("RecordWithGeneric3", "")
			.declare_key("int_array3", Instance(generic_rec, param_vec).close(), "key desc")
			.close();

		generic_rec.finish(FinishStatus::generic_);
		EXPECT_THROW_WHAT( { rec.finish(); }, ExcWrongDefault,
				"do not match type: 'array_of_Integer'" );
	}

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
    EXPECT_EQ("TYPE", it->key_);
    ++it;
    EXPECT_EQ("file", it->key_);
    EXPECT_EQ("File for output stream.", it->description_);
    EXPECT_EQ(typeid(FileName) , typeid(*(it->type_)));
    EXPECT_EQ(FilePath::output_file, static_cast<const FileName *>( &(*it->type_) )->get_file_type() );
    it+=4;
    EXPECT_EQ( "data_description", it->key_);
    EXPECT_EQ( output_record.end(), it+1 );
    // method size()
    EXPECT_EQ(6, output_record.size());
    //method key_index
    EXPECT_EQ(1,output_record.key_index("file"));
    EXPECT_EQ(3,output_record.key_index("compression"));
    EXPECT_THROW({output_record.key_index("x_file");}, Record::ExcRecordKeyNotFound );
    // method key_iterator
    EXPECT_EQ(3,output_record.key_iterator("compression")->key_index);


}

TEST(InputTypeRecord, check_key_validity) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Record output_record("OutputRecord",
            "Information about one file for field data.");
    EXPECT_THROW_WHAT( {output_record.declare_key("a b",Bool(),"desc."); },
    		feal::Exc_assert, "Invalid key identifier"
            );
    EXPECT_THROW_WHAT( {output_record.declare_key("AB",Bool(),"desc."); },
    		feal::Exc_assert, "Invalid key identifier"
            );
    EXPECT_THROW_WHAT( {output_record.declare_key("%$a",Bool(),"desc."); },
    		feal::Exc_assert, "Invalid key identifier"
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

    EXPECT_EQ(4, composite.size());
    EXPECT_EQ("a from rec1", composite.key_iterator("a")->description_);
    EXPECT_EQ("b from composite", composite.key_iterator("b")->description_);
    EXPECT_EQ("c from rec2", composite.key_iterator("c")->description_);
}
