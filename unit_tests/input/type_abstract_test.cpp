/*
 * type_abstract_test.cpp
 *
 *  Created on: May 4, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>

#include <input/type_record.hh>



/**
 * Test Abstract Record.
 */

TEST(InputTypeAbstractRecord, inheritance) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

	Record copy_rec = Record("Copy","")
       	.declare_key("mesh", String(), Default("input.msh"), "Comp. mesh.")
       	.declare_key("a_val", String(), Default::obligatory(), "")
		.close();

    AbstractRecord a_rec = AbstractRecord("EqBase","Base of equation records.");
    AbstractRecord &a_ref = a_rec.allow_auto_conversion("EqDarcy").close();
    EXPECT_EQ( a_rec, a_ref);

    // test derived type
    Record b_rec = Record("EqDarcy","")
    	.derive_from(a_rec)
		.copy_keys(copy_rec)
    	.declare_key("b_val", Integer(), Default("10"), "")
    	.allow_auto_conversion("a_val")
		.close();

    Record c_rec = Record("EqTransp","")
    	.derive_from(a_rec)
		.copy_keys(copy_rec)
    	.declare_key("c_val", Integer(), "")
    	.declare_key("a_val", Double(),"")
		.close();

    a_rec.add_child(b_rec);
    a_rec.add_child(c_rec);
    c_rec.finish();
    b_rec.finish();
    a_rec.finish();

    // auto conversion - default value for TYPE
    EXPECT_EQ("EqDarcy", a_rec.get_selection_default().value() );
    // no more allow_auto_conversion for a_rec
    EXPECT_THROW_WHAT( { a_rec.allow_auto_conversion("EqTransp");}, ExcXprintfMsg, "Can not specify default value for TYPE key as the AbstractRecord 'EqBase' is closed.");

    a_rec.finish();
    EXPECT_EQ( b_rec,  * a_rec.get_default_descendant() );

    // test default value for an auto convertible abstract record key
    Record xx_rec = Record("XX", "")
    		.declare_key("ar_key", a_rec, Default("ahoj"), "")
			.close();
    xx_rec.finish();

    // check correct stat of a_rec
    EXPECT_TRUE( a_rec.is_finished() );
    //EXPECT_EQ(0, a_rec.key_index("TYPE"));
    //EXPECT_EQ(Selection("EqBase_TYPE_selection"), *(a_rec.key_iterator("TYPE")->type_ ));
    //EXPECT_EQ(a_rec.key_iterator("TYPE")->type_, b_rec.key_iterator("TYPE")->type_);

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
    AbstractRecord x = AbstractRecord("AR","")
    	.allow_auto_conversion("BR")
		.close();
    EXPECT_THROW_WHAT({ x.finish(); }, ExcXprintfMsg, "Default value 'BR' for TYPE key do not match any descendant of AbstractRecord 'AR'.");

}





/**
 * Test AdHocAbstractRecord.
 */
/*namespace IT=Input::Type;

class AdHocDataTest : public testing::Test {
public:
	static const IT::Record & get_rec();
	static const IT::Record & get_in_rec1();
	static const IT::Record & get_in_rec2();
	static const IT::AbstractRecord & get_ancestor();
	static const IT::AdHocAbstractRecord & get_adhoc_1();
	static const IT::AdHocAbstractRecord & get_adhoc_2();

protected:
    virtual void SetUp() {
    }
    virtual void TearDown() {
    };
};


const IT::Record & AdHocDataTest::get_in_rec1() {
	return IT::Record("Record 1","")
		.declare_key("val_1", IT::Integer(0), "value 1" )
		.close();
}

const IT::AdHocAbstractRecord & AdHocDataTest::get_adhoc_1() {
	return IT::AdHocAbstractRecord(get_ancestor())
	.add_child(AdHocDataTest::get_in_rec1())
	.add_child(AdHocDataTest::get_in_rec2());
}

const IT::Record & AdHocDataTest::get_rec() {
	return IT::Record("Problem","Base record")
		.declare_key("adhoc_1", AdHocDataTest::get_adhoc_1(), "" )
		.declare_key("adhoc_2", AdHocDataTest::get_adhoc_2(), "" )
		.close();
}

const IT::AdHocAbstractRecord & AdHocDataTest::get_adhoc_2() {
	return IT::AdHocAbstractRecord(get_ancestor())
		.add_child(AdHocDataTest::get_in_rec1())
		.add_child(AdHocDataTest::get_in_rec2());
}

const IT::AbstractRecord & AdHocDataTest::get_ancestor() {
	return IT::AbstractRecord("Ancestor","Base of equation records.").close();
}

const IT::Record & AdHocDataTest::get_in_rec2() {
	return IT::Record("Record 2","")
		.declare_key("val_2", IT::Integer(0), "value 2" )
		.close();
}


TEST(InputTypeAdHocAbstractRecord, inheritance) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";
	AdHocDataTest::get_in_rec1().finish();
	AdHocDataTest::get_in_rec2().finish();
	AdHocDataTest::get_adhoc_1().finish();
	AdHocDataTest::get_adhoc_2().finish();
	AdHocDataTest::get_rec().finish();

	EXPECT_EQ( 1, AdHocDataTest::get_in_rec1().size());
	EXPECT_EQ( 1, AdHocDataTest::get_in_rec2().size());
	EXPECT_EQ( 2, AdHocDataTest::get_adhoc_1().child_size());
	EXPECT_EQ( 2, AdHocDataTest::get_adhoc_2().child_size());
	EXPECT_EQ( 2, AdHocDataTest::get_rec().size());
}*/
