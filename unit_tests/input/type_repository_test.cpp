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
    static const IT::Record get_main_record();
    static const IT::Record get_stat_rec();
    static const IT::Selection get_stat_sel();
    static const IT::AbstractRecord get_stat_a_rec();
};

// test that declare_key method correctly deal with both static (not yet constructed) and local
// types
// The same we test for all types of Array
const IT::Record InputTypeCollection::get_main_record() {
	IT::Record type = IT::Record("MainRecord","")
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
        .declare_key("arr_stat_a_rec", IT::Array(InputTypeCollection::get_stat_a_rec()), "" );

	return type.close();
}




const IT::Record InputTypeCollection::get_stat_rec() {
	return IT::Record("stat_rec","").close();
}

const IT::Selection InputTypeCollection::get_stat_sel() {
	return IT::Selection("stat_sel").close();
}


const IT::AbstractRecord InputTypeCollection::get_stat_a_rec() {
	IT::AbstractRecord a_rec = IT::AbstractRecord("stat_a_rec","");
	a_rec.close();
	return a_rec.close();
}




TEST(TypeRepository, all) {
    IT::Record rec = InputTypeCollection::get_main_record();
    //std::cout << OutputText( &rec ) << std::endl;

    EXPECT_EQ( InputTypeCollection::get_stat_rec(), InputTypeCollection::get_stat_rec() );
    //std::cout << &(InputTypeCollection::get_stat_rec()) << endl;
    //std::cout << &(InputTypeCollection::get_stat_rec()) << endl;

}
