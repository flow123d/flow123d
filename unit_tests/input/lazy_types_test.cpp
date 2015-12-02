/*
 * lazy_types_test.cpp
 *
 *  Created on: Jan 12, 2013
 *      Author: jb
 */



#include <flow_gtest.hh>

#include <fstream>
#include <input/input_type.hh>
#include <input/type_output.hh>

using namespace Input::Type;
namespace IT=Input::Type;

IT::Record get_loc_rec() {
    return IT::Record("loc_rec","")
            .declare_key("key",Integer(),"");
}

IT::Selection get_loc_sel() {
    return IT::Selection("loc_sel")
            .add_value(13,"happy");
}

// test local record derive from local abstract record
IT::Record loc_a_rec_descendant(IT::Abstract &type) {
    return Record("loc_a_rec_descendant","")
            .derive_from(type)
            .declare_key("str", String(),"");
}

IT::Abstract get_loc_a_rec() {
    IT::Abstract type = IT::Abstract("loc_a_rec","")
            .declare_key("key",Integer(),"");

    // local instance of Abstract has to call creation of its descendants
    // in order to make them all derived from very same instance of Abstract
    // different calls to get_loc_a_rec() returns completely distinguish types even if they
    // has same name and structure. In particular no static Record can derive from local Abstract !!

    loc_a_rec_descendant(type);
    return type;
}


struct InputTypeCollection {
    static IT::Record main_record;
    static IT::Record stat_rec;
    static IT::Selection stat_sel;
    static IT::Abstract stat_a_rec;
    static IT::Record stat_desc_stat_a_rec;
};

// test that declare_key method correctly deal with both static (not yet constructed) and local
// types
// The same we test for all types of Array
IT::Record InputTypeCollection::main_record=
        IT::Record("MainRecord","")
        .declare_key("stat_rec", stat_rec, "")
        .declare_key("loc_rec", get_loc_rec(), "")
        .declare_key("loc_rec_2", get_loc_rec(), "")
        .declare_key("stat_sel", stat_sel, "")
        .declare_key("loc_sel", get_loc_sel(), "")
        .declare_key("stat_a_rec", stat_a_rec, "" )
        .declare_key("loc_a_rec", get_loc_a_rec(), "" )

        .declare_key("arr_stat_rec", IT::Array(stat_rec), "")
        .declare_key("arr_loc_rec", IT::Array(get_loc_rec()), "")
        .declare_key("arr_loc_rec_2", IT::Array(get_loc_rec()), "")
        .declare_key("arr_stat_sel", IT::Array(stat_sel), "")
        .declare_key("arr_loc_sel", IT::Array(get_loc_sel()), "")
        .declare_key("arr_stat_a_rec", IT::Array(stat_a_rec), "" )
        .declare_key("arr_loc_a_rec", IT::Array(get_loc_a_rec()), "" );




IT::Record InputTypeCollection::stat_rec=IT::Record("stat_rec","");

IT::Selection InputTypeCollection::stat_sel=IT::Selection("stat_sel");

// test static record derive from static abstract record
IT::Record InputTypeCollection::stat_desc_stat_a_rec=
        IT::Record("stat_desc_stat_a_rec","")
        .derive_from(stat_a_rec);


IT::Abstract InputTypeCollection::stat_a_rec=IT::Abstract("stat_a_rec","");




TEST(LazyTypes, all) {
    TypeBase::lazy_finish();

    std::cout << OutputText(&InputTypeCollection::main_record) << std::endl;

}
