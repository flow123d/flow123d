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
//::testing::FLAGS_gtest_death_test_style = "threadsafe";


   // make auxiliary record and test declare_key for
   // - various Scalar types (excluding Selection): Integer, Bool, Double, String, FilePath
   // - various decalre_key templates: with/without default, shared_ptr/ reference
   // - various default values
   static Record rec("SomeRecord", "desc.");

   rec.declare_key("file", FileName::output(), Default::optional(), "desc.");

   Integer digits_type(0, 8);
   rec.declare_key("digits",digits_type, Default("8"), "desc.");

   rec.declare_key("compression", Bool(),"desc.");

   Double time_type(0.0);
   rec.declare_key("start_time", time_type,"desc.");

   rec.declare_key("data_description", String(), Default::optional(),"");


   // errors during declaration
   Record rec_empty;

#ifdef DEBUG_ASSERTS
   EXPECT_THROW_WHAT( {rec_empty.declare_key("xx", Integer(), "");}, ExcXprintfMsg, "Internal Error"); //".*into closed record 'EmptyRecord'."
#endif

   Record rec_fin("xx","");
   rec_fin.close();
   EXPECT_THROW_WHAT( {rec_fin.declare_key("xx", String(),"");}, ExcXprintfMsg, "Internal Error"); //"Can not add .* into closed record"


//   This no more fails: Declaration of incomplete (unfinished) keys is possible.
   Record rec_unfin("yy","");
   rec.declare_key("yy", rec_unfin, "");

   EXPECT_THROW_WHAT( { rec.declare_key("data_description", String(),"");}, ExcXprintfMsg, "Error"); //"Re-declaration of the key:"

   EXPECT_THROW_WHAT( { rec.declare_key("wrong_double", Double(), Default("1.23 4"),""); }, ExcWrongDefault,
           "Default value .* do not match type: 'Double';");

   enum Colors {
       white, black, red
   };

   Selection sel("Color selection");
   sel.add_value(black, "black");
   sel.add_value(red, "red");
   sel.close();

   rec.declare_key("plot_color", sel, "Color to plot the fields in file.");

   // test correct finishing.
#ifdef DEBUG_ASSERTS
   EXPECT_THROW_WHAT( {rec.size();}, ExcAssertMsg , "Asking for information of unfinished Record type");
#endif

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


   static Record array_record("RecordOfArrays", "desc.");
   static Record array_record2("RecordOfArrays2", "desc.");

   // array type passed through shared_ptr
    Array array_of_int(Integer(0), 5, 100 );
    array_record.declare_key("array_of_5_ints", array_of_int,"Some bizare array.");

    // array type passed by reference
    array_record.declare_key("array_of_str", Array( String() ),"Desc. of array");
    array_record.declare_key("array_of_str_1", Array( String() ), "Desc. of array");


    // allow default values for an array
    array_record.declare_key("array_with_default", Array( Double() ), Default("3.2"), "");
    EXPECT_THROW_WHAT( {array_record.declare_key("some_key", Array( Integer() ), Default("ahoj"), ""); array_record.close(); }, ExcWrongDefault,
                  "Default value 'ahoj' do not match type: 'Integer';"
                 );
    EXPECT_THROW_WHAT( {array_record2.declare_key("some_key", Array( Double(), 2 ), Default("3.2"), ""); array_record2.close(); }, ExcWrongDefault,
                  "Default value '3.2' do not match type: 'array_of_Double';"
                 );
}

TEST(InputTypeRecord, allow_convertible) {
using namespace Input::Type;
//::testing::FLAGS_gtest_death_test_style = "threadsafe";

    {
    static Record sub_rec( "SubRecord", "");
    sub_rec.declare_key("bool_key", Bool(), Default("false"), "");
    sub_rec.declare_key("int_key", Integer(),  "");
    sub_rec.allow_auto_conversion("int_key");
    sub_rec.finish();

    EXPECT_EQ(1, sub_rec.auto_conversion_key_iter()->key_index );
    }

    {
    static Record sub_rec( "SubRecord", "");
    sub_rec.declare_key("bool_key", Bool(), "");
    sub_rec.declare_key("int_key", Integer(),  "");
    sub_rec.allow_auto_conversion("int_key");
    EXPECT_THROW_WHAT( {sub_rec.finish();}, ExcXprintfMsg, "Internal Error"); //"Finishing Record auto convertible from the key 'int_key', but other key: 'bool_key' has no default value."
    }

    {
    static Record sub_rec( "SubRecord", "");
    EXPECT_TRUE( sub_rec.finish() );
    EXPECT_EQ( sub_rec.end(), sub_rec.auto_conversion_key_iter() );
    }

}


TEST(InputTypeRecord, declare_key_record) {
using namespace Input::Type;
//::testing::FLAGS_gtest_death_test_style = "threadsafe";


    Record record_record("RecordOfRecords", "");
    Record record_record2("RecordOfRecords2", "");

    // Test that Record has to be passed as shared_ptr
    //ASSERT_DEATH( {record_record->declare_key("sub_rec_1", Record("subrec_type", "desc") , "desc"); },
    //              "Complex type .* shared_ptr."
    //              );

    Record other_record("OtherRecord","desc");
    other_record.close();

    static Record sub_rec( "SubRecord", "");
    sub_rec.declare_key("bool_key", Bool(), Default("false"), "");
    sub_rec.declare_key("int_key", Integer(),  "");
    sub_rec.allow_auto_conversion("int_key");
    sub_rec.close();

    record_record.declare_key("sub_rec_1", other_record, "key desc");
    EXPECT_THROW_WHAT( { record_record.declare_key("sub_rec_2", other_record, Default("2.3"), "key desc"); record_record.close(); }, ExcWrongDefault,
            "Default value '2.3' do not match type: 'OtherRecord';" );

    DBGMSG("here\n");
    record_record2.declare_key("sub_rec_dflt", sub_rec, Default("123"), "");
    DBGMSG("here\n");
    EXPECT_THROW_WHAT( { record_record2.declare_key("sub_rec_dflt2", sub_rec, Default("2.3"), ""); record_record2.close(); } , ExcWrongDefault,
            "Default value '2.3' do not match type: 'Integer';" );

    // recursion  -  forbidden
    //record_record->declare_key("sub_rec_2", record_record, "desc.");

}

TEST(InputTypeRecord, iterating) {
    using namespace Input::Type;

    Record output_record("OutputRecord",
            "Information about one file for field data.");
    {
        output_record.declare_key("file", FileName::output(), Default::optional(),
                "File for output stream.");

        Integer digits_type((int)0, (int)8);
        output_record.declare_key("digits",digits_type, Default("8"),
                "Number of digits used for output double values into text output files.");
        output_record.declare_key("compression", Bool(),
                "Whether to use compression of output file.");

        Double time_type(0.0);
        output_record.declare_key("start_time", time_type,
                "Simulation time of first output.");
        output_record.declare_key("data_description", String(),
                "");
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
//::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Record output_record("OutputRecord",
            "Information about one file for field data.");
    EXPECT_THROW_WHAT( {output_record.declare_key("a b",Bool(),"desc."); },
    		ExcXprintfMsg, "Internal Error" //"Invalid key identifier"
            );
    EXPECT_THROW_WHAT( {output_record.declare_key("AB",Bool(),"desc."); },
    		ExcXprintfMsg, "Internal Error" //"Invalid key identifier"
            );
    EXPECT_THROW_WHAT( {output_record.declare_key("%$a",Bool(),"desc."); },
    		ExcXprintfMsg, "Internal Error" //"Invalid key identifier"
            );

}

TEST(InputTypeRecord, RecordCopy) {
using namespace Input::Type;
//::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Record output_record("OutputRecord", "");
    output_record.declare_key("file", FileName::output(), "");


    Record copy_rec = output_record;

    Integer digits_type((int)0, (int)8);
    copy_rec.declare_key("digits",digits_type, "");
    copy_rec.finish();

    EXPECT_EQ( true, output_record.is_finished());
    EXPECT_EQ( true, output_record.has_key("digits") );
}

/**
 * Test Abstract Record.
 */

TEST(InputTypeAbstractRecord, inheritance) {
using namespace Input::Type;
//::testing::FLAGS_gtest_death_test_style = "threadsafe";

    AbstractRecord a_rec("EqBase","Base of equation records.");
    a_rec.declare_key("mesh", String(), Default("input.msh"), "Comp. mesh.");
    a_rec.declare_key("a_val", String(), Default::obligatory(), "");
    AbstractRecord &a_ref = a_rec.allow_auto_conversion("EqDarcy");
    EXPECT_EQ( a_rec, a_ref);
    a_rec.finish();

    // test derived type
    Record b_rec("EqDarcy","");
    b_rec.derive_from(a_rec);
    b_rec.declare_key("b_val", Integer(), Default("10"), "");
    b_rec.allow_auto_conversion("a_val");

    Record c_rec("EqTransp","");
    c_rec.derive_from(a_rec);
    c_rec.declare_key("c_val", Integer(), "");
    c_rec.declare_key("a_val", Double(),"");

    c_rec.finish();
    b_rec.finish();

    // auto conversion - default value for TYPE
    EXPECT_EQ("EqDarcy", a_rec.key_iterator("TYPE")->default_.value() );
    // no more allow_auto_conversion for a_rec
    EXPECT_THROW_WHAT( { a_rec.allow_auto_conversion("EqTransp");}, ExcXprintfMsg, "Internal Error"); //"Can not specify default value for TYPE key as the AbstractRecord 'EqBase' is closed."

    a_rec.no_more_descendants();
    EXPECT_EQ( b_rec,  * a_rec.get_default_descendant() );

    // test default value for an auto convertible abstract record key
    Record xx_rec("XX", "");
    xx_rec.declare_key("ar_key", a_rec, Default("ahoj"), "");
    xx_rec.finish();

    // check correct stat of a_rec
    EXPECT_TRUE( a_rec.is_finished() );
    EXPECT_EQ(0, a_rec.key_index("TYPE"));
    EXPECT_EQ(Selection("EqBase_TYPE_selection"), *(a_rec.key_iterator("TYPE")->type_ ));
    EXPECT_EQ(a_rec.key_iterator("TYPE")->type_, b_rec.key_iterator("TYPE")->type_);

    // TYPE should be derived as optional
    //EXPECT_TRUE( b_rec.key_iterator("TYPE")->default_.is_optional());
    //EXPECT_TRUE( c_rec.key_iterator("TYPE")->default_.is_optional());

    // inherited keys
    EXPECT_TRUE( b_rec.has_key("mesh") );
    EXPECT_TRUE( c_rec.has_key("mesh") );
    // overwritten key
    EXPECT_EQ( Double(), *(c_rec.key_iterator("a_val")->type_));

    //get descendant
    EXPECT_EQ( b_rec, a_rec.get_descendant("EqDarcy"));
    EXPECT_EQ( c_rec, a_rec.get_descendant("EqTransp"));


    // check of correct auto conversion value
    AbstractRecord  x("AR","");
    x.allow_auto_conversion("BR");
    EXPECT_THROW_WHAT({ x.no_more_descendants(); }, ExcXprintfMsg, "Internal Error"); //"Default value 'BR' for TYPE key do not match any descendant of AbstractRecord 'AR'."

}


/**
 * Test AdHocAbstractRecord.
 */
namespace IT=Input::Type;

class AdHocDataTest : public testing::Test {
public:
	static IT::Record rec;
	static IT::Record in_rec1;
	static IT::Record in_rec2;
	static IT::AbstractRecord ancestor;
	static IT::AdHocAbstractRecord adhoc_1;
	static IT::AdHocAbstractRecord adhoc_2;

protected:
    virtual void SetUp() {
    }
    virtual void TearDown() {
    };
};


IT::Record AdHocDataTest::in_rec1 = IT::Record("Record 1","")
	.declare_key("val_1", IT::Integer(0), "value 1" )
	.close();

IT::AdHocAbstractRecord AdHocDataTest::adhoc_1 = IT::AdHocAbstractRecord(ancestor)
	.add_child(AdHocDataTest::in_rec1)
	.add_child(AdHocDataTest::in_rec2);

IT::Record AdHocDataTest::rec = IT::Record("Problem","Base record")
	.declare_key("adhoc_1", AdHocDataTest::adhoc_1, "" )
	.declare_key("adhoc_2", AdHocDataTest::adhoc_2, "" );

IT::AdHocAbstractRecord AdHocDataTest::adhoc_2 = IT::AdHocAbstractRecord(ancestor)
	.add_child(AdHocDataTest::in_rec1)
	.add_child(AdHocDataTest::in_rec2);

IT::AbstractRecord AdHocDataTest::ancestor = IT::AbstractRecord("Ancestor","Base of equation records.");

IT::Record AdHocDataTest::in_rec2 = IT::Record("Record 2","")
	.declare_key("val_2", IT::Integer(0), "value 2" )
	.close();


TEST(InputTypeAdHocAbstractRecord, inheritance) {
using namespace Input::Type;
//::testing::FLAGS_gtest_death_test_style = "threadsafe";
	AdHocDataTest::in_rec1.finish();
	AdHocDataTest::in_rec2.finish();
	AdHocDataTest::adhoc_1.finish();
	AdHocDataTest::adhoc_2.finish();
	AdHocDataTest::rec.finish();

	EXPECT_EQ( 1, AdHocDataTest::in_rec1.size());
	EXPECT_EQ( 1, AdHocDataTest::in_rec2.size());
	EXPECT_EQ( 2, AdHocDataTest::adhoc_1.child_size());
	EXPECT_EQ( 2, AdHocDataTest::adhoc_2.child_size());
	EXPECT_EQ( 2, AdHocDataTest::rec.size());
}
