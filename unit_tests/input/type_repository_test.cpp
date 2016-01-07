/*
 * type_repository_test.cpp
 *
 *  Created on: Jan 12, 2013
 *      Author: jb
 */



#include <flow_gtest.hh>

#include <fstream>
#include <input/input_type.hh>
#include <input/type_output.hh>
#include <input/type_repository.hh>

using namespace Input::Type;
namespace IT=Input::Type;

IT::Record get_loc_rec() {
    return IT::Record("loc_rec","")
            .declare_key("key",Integer(),"").close();
}

IT::Selection get_loc_sel() {
    return IT::Selection("loc_sel")
            .add_value(13,"happy").close();
}


struct InputTypeCollection {
    static const IT::Record & get_main_record();
    static const IT::Record & get_stat_rec();
    static const IT::Selection & get_stat_sel();
    static const IT::Abstract & get_stat_a_rec();
};

// test that declare_key method correctly deal with both static (not yet constructed) and local
// types
// The same we test for all types of Array
const IT::Record & InputTypeCollection::get_main_record() {
	return IT::Record("MainRecord","")
        .declare_key("stat_rec", InputTypeCollection::get_stat_rec(), "")
        .declare_key("loc_rec", get_loc_rec(), "")
        .declare_key("loc_rec_2", get_loc_rec(), "")
        .declare_key("stat_sel", InputTypeCollection::get_stat_sel(), "")
        .declare_key("loc_sel", get_loc_sel(), "")
        .declare_key("stat_a_rec", InputTypeCollection::get_stat_a_rec(), "" )

        .declare_key("arr_stat_rec", IT::Array(InputTypeCollection::get_stat_rec()), "")
        .declare_key("arr_loc_rec", IT::Array(get_loc_rec()), "")
        .declare_key("arr_loc_rec_2", IT::Array(get_loc_rec()), "")
        .declare_key("arr_stat_sel", IT::Array(InputTypeCollection::get_stat_sel()), "")
        .declare_key("arr_loc_sel", IT::Array(get_loc_sel()), "")
        .declare_key("arr_stat_a_rec", IT::Array(InputTypeCollection::get_stat_a_rec()), "" )
		.close();
}




const IT::Record & InputTypeCollection::get_stat_rec() {
	return IT::Record("stat_rec","").close();
}

const IT::Selection & InputTypeCollection::get_stat_sel() {
	return IT::Selection("stat_sel").close();
}


const IT::Abstract & InputTypeCollection::get_stat_a_rec() {
	return IT::Abstract("stat_a_rec","").close();
}




TEST(TypeRepository, all) {
    IT::Record rec = InputTypeCollection::get_main_record();

    EXPECT_EQ( InputTypeCollection::get_stat_rec(), InputTypeCollection::get_stat_rec() );

}





TEST(TypeRepository, hash) {
	TypeBase::TypeHash hash;

	// test of Record
	IT::Record rec = IT::Record("record", "Some record.")
		.declare_key("int_key", Integer(),"")
		.declare_key("double_key", Double(),"");

	hash = rec.content_hash();
	EXPECT_EQ( hash, rec.close().content_hash() ); // close() method get record from TypeRepository
	EXPECT_EQ( IT::Record("stat_rec","").close().content_hash(), InputTypeCollection::get_stat_rec().content_hash() );

	// test of Abstract
	IT::Abstract a_rec = IT::Abstract("abstract_record", "Some abstract record.");

	hash = a_rec.content_hash();
	EXPECT_EQ( hash, a_rec.close().content_hash() );
	EXPECT_EQ( IT::Abstract("stat_a_rec","").close().content_hash(), InputTypeCollection::get_stat_a_rec().content_hash() );

	// test of Selection
	IT::Selection sel = IT::Selection("selection", "Some selection.")
		.add_value(0, "off")
		.add_value(1, "on");

	hash = sel.content_hash();
	EXPECT_EQ( hash, sel.close().content_hash() );
	EXPECT_EQ( IT::Selection("stat_sel").close().content_hash(), InputTypeCollection::get_stat_sel().content_hash() );

}

