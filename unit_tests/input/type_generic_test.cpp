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

	static Record rec = Record("generic_rec", "desc.")
							.root_of_generic_subtree()
							.declare_key("param1", Parameter("param1"), "desc.")
							.declare_key("param2", Parameter("param2"), "desc.")
							.declare_key("start_time", Double(), "desc.")
							.declare_key("name", String(), "desc.")
							.close();

	return Instance(rec, param_vec)
			.close();
}

static const Instance & get_generic_array(const Selection *sel) {
	static Array arr = Array(Parameter("param"), 0, 100).root_of_generic_subtree();

	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param", boost::make_shared<Selection>(*sel)) );

	return Instance( arr, param_vec )
			.close();
}

static AbstractRecord & get_abstract_type() {
	return AbstractRecord("SomeAbstract", "Some Abstract.").root_of_generic_subtree().close();
}

static const Instance & get_generic_abstract(const Selection *sel) {
	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param", boost::make_shared<Selection>(*sel)) );

	return Instance( get_abstract_type(), param_vec )
			.close();
}

static const Record & get_inner_record() {
	return Record("inner_rec", "desc.")
			.declare_key("param2", Parameter("param2"), "desc.")
			.declare_key("name", String(), "desc.")
			.close();
}

/**
 * Returned part of IST in format:
 *
 * Record
 *  - Record
 *     - Parameter
 *     - String
 *  - Parameter
 *  - Integer
 */
static const Instance & get_record_with_record(const Selection *sel, const boost::shared_ptr<TypeBase> type) {
	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param1", boost::make_shared<Selection>(*sel)) );
	param_vec.push_back( std::make_pair("param2", type) );

	static Record rec = Record("record", "desc.")
					.root_of_generic_subtree()
					.declare_key("inner_record", get_inner_record(), "desc.")
					.declare_key("param1", Parameter("param1"), "desc.")
					.declare_key("start_time", Double(), "desc.")
					.close();

	return Instance(rec, param_vec)
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


TEST(GenericType, record_with_record) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";
	static Record problem = Record("problem", "desc.")
			.declare_key("primary", get_record_with_record(&get_colors_selection(), boost::make_shared<Integer>(0,100)), "Primary problem.")
			.declare_key("secondary", get_record_with_record(&get_shapes_selection(), boost::make_shared<Double>()), "Secondary problem.")
			.declare_key("bool", Bool(), "Some bool key.")
			.close();

	TypeBase::lazy_finish();

	Record::KeyIter key_it = problem.begin();
	{
		EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Record) );
		const Record *in_rec = static_cast<const Record *>(key_it->type_.get());
		in_rec->write_attributes(cout); cout << endl;
		Record::KeyIter in_rec_it = in_rec->begin();
		EXPECT_EQ( typeid( *(in_rec_it->type_.get()) ), typeid(Record) );
		const Record *in_in_rec = static_cast<const Record *>(in_rec_it->type_.get());
		EXPECT_EQ( typeid( *(in_in_rec->begin()->type_.get()) ), typeid(Integer) );
	}
	++key_it;
	{
		EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Record) );
		const Record *in_rec = static_cast<const Record *>(key_it->type_.get());
		in_rec->write_attributes(cout); cout << endl;
		Record::KeyIter in_rec_it = in_rec->begin();
		EXPECT_EQ( typeid( *(in_rec_it->type_.get()) ), typeid(Record) );
		const Record *in_in_rec = static_cast<const Record *>(in_rec_it->type_.get());
		EXPECT_EQ( typeid( *(in_in_rec->begin()->type_.get()) ), typeid(Double) );
	}
}


TEST(GenericType, parameter_in_deep) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param", boost::make_shared<Double>()) );

	static Record with_array = Record("inner_rec", "")
			.root_of_generic_subtree()
			.declare_key("array", Array( Array( Parameter("param") ) ), "desc.")
			.declare_key("some_double", Double(), "Double key")
			.close();

	static Instance inst = Instance(with_array, param_vec).close();

	static Record record = Record("parametrized_record", "")
			.declare_key("generic_rec", inst, "desc.")
			.declare_key("some_int", Integer(), "Int key")
			.close();

	TypeBase::lazy_finish();
}


TEST(GenericType, array_of_instances) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param", boost::make_shared<Double>()) );

	static Parameter param = Parameter("param");

	static Instance inst = Instance(param, param_vec).close();

	static Array array = Array( inst );

	array.finish();
	TypeBase::lazy_finish();

	EXPECT_EQ( typeid( array.get_sub_type() ), typeid(Double) );
}


TEST(GenericType, instance_in_instance) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	std::vector<TypeBase::ParameterPair> param_vec1, param_vec2;
	param_vec1.push_back( std::make_pair("param", boost::make_shared<Double>()) );
	param_vec2.push_back( std::make_pair("param", boost::make_shared<String>()) );

	static Parameter param = Parameter("param");
	static Instance inst_in = Instance(param, param_vec2).close();

	static Record in_rec = Record("Inner record", "")
			.root_of_generic_subtree()
			.declare_key("key1", inst_in, "Inner instance.")
			.declare_key("param", Parameter("param"), "Parameterized key")
			.close();

	static Instance inst = Instance(in_rec, param_vec1).close();

	static Record root_rec = Record("Root record", "")
			.declare_key("rec", inst, "First instance.")
			.declare_key("some_int", Integer(), "Int key")
			.close();

	TypeBase::lazy_finish();

	Record::KeyIter key_it = root_rec.begin();
	EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Record) );
	const Record *inner_rec = static_cast<const Record *>(key_it->type_.get());
	EXPECT_EQ(inner_rec->size(), 2);
	Record::KeyIter key_it_in = inner_rec->begin();
	EXPECT_EQ( typeid( *(key_it_in->type_.get()) ), typeid(String) );
	++key_it_in;
	EXPECT_EQ( typeid( *(key_it_in->type_.get()) ), typeid(Double) );
}


TEST(GenericType, parameter_not_replaced) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	std::vector<TypeBase::ParameterPair> param_vec;

	static Record inner = Record("inner_rec", "")
			.root_of_generic_subtree()
			.declare_key("param", Parameter("param"), "desc.")
			.declare_key("some_double", Double(), "Double key")
			.close();

	static Instance inst = Instance(inner, param_vec)
								.close();

	static Record record = Record("parametrized_record", "")
			.declare_key("generic_rec", inst, "desc.")
			.declare_key("some_double", Double(), "Double key")
			.close();

	EXPECT_THROW_WHAT( { TypeBase::lazy_finish(); }, ExcParamaterNotSubsituted, "for input type parameter 'param'" );
}


TEST(GenericType, parameter_not_used) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param", boost::make_shared<Integer>()) );
	param_vec.push_back( std::make_pair("param2", boost::make_shared<Double>()) );

	static Record inner = Record("inner_rec", "")
			.root_of_generic_subtree()
			.declare_key("param", Parameter("param"), "desc.")
			.declare_key("some_double", Double(), "Double key")
			.close();

	static Instance inst = Instance(inner, param_vec)
								.close();

	static Record record = Record("parametrized_record", "")
			.declare_key("generic_rec", inst, "desc.")
			.declare_key("some_double", Double(), "Double key")
			.close();

	EXPECT_ASSERT_DEATH( { TypeBase::lazy_finish(); }, "must be used");
}


TEST(GenericType, parameter_during_lazy_finish) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	static Record record = Record("parametrized_record", "")
			.declare_key("param", Parameter("param"), "desc.")
			.declare_key("some_double", Double(), "Double key")
			.close();

	EXPECT_ASSERT_DEATH( { TypeBase::lazy_finish(); }, "Finish of non-generic Parameter 'param'");
}
