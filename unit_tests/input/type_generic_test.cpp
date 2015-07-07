/*
 * input_generic_test.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *
 *
 */


#include <flow_gtest.hh>

#include <input/type_generic.hh>
#include <input/type_base.hh>
#include <input/type_record.hh>


using namespace Input::Type;



static const Selection & get_colors_selection() {
	return Selection("colors")
			.add_value(0, "red")
			.add_value(1, "green")
			.add_value(2, "blue")
			.close();
}

static const Selection & get_shapes_selection() {
	return Selection("shapes")
			.add_value(0, "rectangle")
			.add_value(1, "circle")
			.add_value(2, "triangle")
			.close();
}

static const Instance & get_generic_record(const Selection *sel, int max_limit) {
	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param1", boost::make_shared<Selection>(*sel)) );
	param_vec.push_back( std::make_pair("param2", boost::make_shared<Integer>(0, max_limit)) );

	return Instance(Record("generic_rec", "desc.")
						.declare_key("param1", Parameter("param1"), "desc.")
						.declare_key("param2", Parameter("param2"), "desc.")
						.declare_key("start_time", Double(), "desc.")
						.declare_key("name", String(), "desc.")
						.close(),
					param_vec)
			.close();
}

static const Instance & get_generic_array(const Selection *sel) {
	static Array arr(Parameter("param"), 0, 100);

	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param", boost::make_shared<Selection>(*sel)) );

	return Instance( arr, param_vec )
			.close();
}

static AbstractRecord & get_abstract_type() {
	return AbstractRecord("SomeAbstract", "Some Abstract.").close();
}

static const Instance & get_generic_abstract(const Selection *sel) {
	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param", boost::make_shared<Selection>(*sel)) );

	return Instance( get_abstract_type(), param_vec )
			.close();
}


TEST(GenericType, generic_record) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";
	static Record problem = Record("problem", "desc.")
			.declare_key("primary", get_generic_record(&get_colors_selection(), 10), "Primary problem.")
			.declare_key("secondary", get_generic_record(&get_shapes_selection(), 1000), "Secondary problem.")
			.declare_key("bool", Bool(), "Some bool key.")
			.close();

	TypeBase::lazy_finish();

	Record::KeyIter key_it = problem.begin();
	EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Record) );
	const Record *in_rec = static_cast<const Record *>(key_it->type_.get());
	EXPECT_EQ( typeid( *(in_rec->begin()->type_.get()) ), typeid(Selection) );

	++key_it;
	EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Record) );
}


TEST(GenericType, generic_array) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";
	static Record arr_rec = Record("rec_of_arr", "desc.")
			.declare_key("array1", get_generic_array(&get_colors_selection()), "Primary problem.")
			.declare_key("array2", get_generic_array(&get_shapes_selection()), "Secondary problem.")
			.close();

	TypeBase::lazy_finish();

	Record::KeyIter key_it = arr_rec.begin();
	EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Array) );
	const Array *array = static_cast<const Array *>(key_it->type_.get());
	EXPECT_EQ( typeid( array->get_sub_type() ), typeid(Selection) );

	++key_it;
	EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Array) );
}


TEST(GenericType, generic_abstract) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";
	static Record desc1 = Record("Descendant1", "")
			.derive_from(get_abstract_type())
			.declare_key("param", Parameter("param"), "desc.")
			.declare_key("some_int", Integer(), "Integer key")
			.close();

	static Record desc2 = Record("Descendant2", "")
			.derive_from(get_abstract_type())
			.declare_key("param", Parameter("param"), "desc.")
			.declare_key("some_double", Double(), "Double key")
			.close();

	static Record rec_with_abstracts = Record("rec_of_arr", "desc.")
			.declare_key("abstract1", get_generic_abstract(&get_colors_selection()), "Primary problem.")
			.declare_key("abstract2", get_generic_abstract(&get_shapes_selection()), "Secondary problem.")
			.declare_key("bool", Bool(), "Some bool key.")
			.close();

	TypeBase::lazy_finish();

	Record::KeyIter key_it = rec_with_abstracts.begin();
	{
		EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(AbstractRecord) );
		const AbstractRecord *abstract = static_cast<const AbstractRecord *>(key_it->type_.get());
		EXPECT_EQ( abstract->child_size(), 2 );
		Record::KeyIter desc_it = abstract->get_descendant(0).begin();
		EXPECT_EQ( typeid( *(desc_it->type_.get()) ), typeid(Selection) );
		++desc_it;
		const Selection *sel = static_cast<const Selection *>(desc_it->type_.get());
		EXPECT_EQ( sel->type_name(), "colors" );
	}

	++key_it;
	{
		EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(AbstractRecord) );
		const AbstractRecord *abstract = static_cast<const AbstractRecord *>(key_it->type_.get());
		EXPECT_EQ( abstract->child_size(), 2 );
		Record::KeyIter desc_it = abstract->get_descendant(0).begin();
		EXPECT_EQ( typeid( *(desc_it->type_.get()) ), typeid(Selection) );
		++desc_it;
		const Selection *sel = static_cast<const Selection *>(desc_it->type_.get());
		EXPECT_EQ( sel->type_name(), "shapes" );
	}
}
