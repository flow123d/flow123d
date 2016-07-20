/*
 * input_generic_test.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *
 *
 */


#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include <input/input_type.hh>


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
	param_vec.push_back( std::make_pair("param1", std::make_shared<Selection>(*sel)) );
	param_vec.push_back( std::make_pair("param2", std::make_shared<Integer>(0, max_limit)) );

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

static const Instance & get_generic_tuple(const std::shared_ptr<TypeBase> type, int max_limit) {
	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param1", type ) );
	param_vec.push_back( std::make_pair("param2", std::make_shared<Integer>(0, max_limit)) );

	static Tuple tuple = Tuple("generic_tuple", "desc.")
							.root_of_generic_subtree()
							.declare_key("param1", Parameter("param1"), Default::obligatory(), "desc.")
							.declare_key("param2", Parameter("param2"), Default::obligatory(), "desc.")
							.declare_key("time", Double(), Default::optional(), "desc.")
							.close();

	return Instance(tuple, param_vec)
			.close();
}

static const Instance & get_generic_array(const Selection *sel) {
	static Array arr = Array(Parameter("param"), 0, 100);

	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param", std::make_shared<Selection>(*sel)) );

	return Instance( arr, param_vec )
			.close();
}

static Abstract & get_abstract_type() {
	return Abstract("SomeAbstract", "Some Abstract.").root_of_generic_subtree().close();
}

static const Instance & get_generic_abstract(const Selection *sel) {
	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param", std::make_shared<Selection>(*sel)) );

	return Instance( get_abstract_type(), param_vec )
			.close();
}

static const Record & get_inner_record() {
	return Record("inner_rec", "desc.")
			.declare_key("param2", Parameter("param2"), "desc.")
			.declare_key("name", String(), "desc.")
			.close();
}

static const Tuple & get_inner_tuple() {
	return Tuple("tuple", "desc.")
			.declare_key("int_key1", Integer(), Default::obligatory(), "desc.")
			.declare_key("int_key2", Integer(), Default("1"), "desc.")
			.declare_key("int_key3", Integer(), Default::optional(), "desc.")
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
 *  - Tuple
 *  - Double
 */
static const Instance & get_record_with_record(const Selection *sel, const std::shared_ptr<TypeBase> type) {
	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param1", std::make_shared<Selection>(*sel)) );
	param_vec.push_back( std::make_pair("param2", type) );

	static Record rec = Record("record", "desc.")
					.root_of_generic_subtree()
					.declare_key("inner_record", get_inner_record(), "desc.")
					.declare_key("param1", Parameter("param1"), "desc.")
					.declare_key("inner_tuple", get_inner_tuple(), "desc.")
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


TEST(GenericType, generic_tuple) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";
	static Record tpl_problem = Record("problem_with_tuple", "desc.")
			.declare_key("primary", get_generic_tuple(std::make_shared<Double>(0.0, 1.0), 10), "Primary problem.")
			.declare_key("secondary", get_generic_tuple(std::make_shared<Double>(1.0), 1000), "Secondary problem.")
			.declare_key("description", String(), "desc.")
			.close();

	TypeBase::lazy_finish();

	Record::KeyIter key_it = tpl_problem.begin();
	{
		EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Tuple) );
		const Tuple *in_tpl = static_cast<const Tuple *>(key_it->type_.get());
		Record::KeyIter in_it = in_tpl->begin();
		EXPECT_EQ( typeid( *(in_it->type_.get()) ), typeid(Double) );
		++in_it;
		EXPECT_EQ( typeid( *(in_it->type_.get()) ), typeid(Integer) );
		++in_it;
		EXPECT_EQ( typeid( *(in_it->type_.get()) ), typeid(Double) );
	}

	++key_it;
	{
		EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Tuple) );
		const Tuple *in_tpl = static_cast<const Tuple *>(key_it->type_.get());
		Record::KeyIter in_it = in_tpl->begin();
		EXPECT_EQ( typeid( *(in_it->type_.get()) ), typeid(Double) );
		++in_it;
		EXPECT_EQ( typeid( *(in_it->type_.get()) ), typeid(Integer) );
		++in_it;
		EXPECT_EQ( typeid( *(in_it->type_.get()) ), typeid(Double) );
	}
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
		EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Abstract) );
		const Abstract *abstract = static_cast<const Abstract *>(key_it->type_.get());
		EXPECT_EQ( abstract->child_size(), 2 );
		Record::KeyIter desc_it = abstract->get_descendant("Descendant1").begin();
		EXPECT_EQ( typeid( *(desc_it->type_.get()) ), typeid(String) );
		++desc_it;
		const Selection *sel = static_cast<const Selection *>(desc_it->type_.get());
		EXPECT_EQ( sel->type_name(), "colors" );
	}

	++key_it;
	{
		EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Abstract) );
		const Abstract *abstract = static_cast<const Abstract *>(key_it->type_.get());
		EXPECT_EQ( abstract->child_size(), 2 );
		Record::KeyIter desc_it = abstract->get_descendant("Descendant1").begin();
		EXPECT_EQ( typeid( *(desc_it->type_.get()) ), typeid(String) );
		++desc_it;
		const Selection *sel = static_cast<const Selection *>(desc_it->type_.get());
		EXPECT_EQ( sel->type_name(), "shapes" );
	}
}


TEST(GenericType, record_with_record) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";
	static Record problem = Record("problem", "desc.")
			.declare_key("primary", get_record_with_record(&get_colors_selection(), std::make_shared<Integer>(0,100)), "Primary problem.")
			.declare_key("secondary", get_record_with_record(&get_shapes_selection(), std::make_shared<Double>()), "Secondary problem.")
			.declare_key("bool", Bool(), "Some bool key.")
			.close();

	TypeBase::lazy_finish();

	Record::KeyIter key_it = problem.begin();
	{
		EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Record) );
		const Record *in_rec = static_cast<const Record *>(key_it->type_.get());
		Record::KeyIter in_rec_it = in_rec->begin();
		EXPECT_EQ( typeid( *(in_rec_it->type_.get()) ), typeid(Record) );
		const Record *in_in_rec = static_cast<const Record *>(in_rec_it->type_.get());
		EXPECT_EQ( typeid( *(in_in_rec->begin()->type_.get()) ), typeid(Integer) );
	}
	++key_it;
	{
		EXPECT_EQ( typeid( *(key_it->type_.get()) ), typeid(Record) );
		const Record *in_rec = static_cast<const Record *>(key_it->type_.get());
		Record::KeyIter in_rec_it = in_rec->begin();
		EXPECT_EQ( typeid( *(in_rec_it->type_.get()) ), typeid(Record) );
		const Record *in_in_rec = static_cast<const Record *>(in_rec_it->type_.get());
		EXPECT_EQ( typeid( *(in_in_rec->begin()->type_.get()) ), typeid(Double) );
	}
}


TEST(GenericType, parameter_in_deep) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param", std::make_shared<Double>()) );

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
	param_vec.push_back( std::make_pair("param", std::make_shared<Double>()) );

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
	param_vec1.push_back( std::make_pair("param", std::make_shared<Double>()) );
	param_vec2.push_back( std::make_pair("param", std::make_shared<String>()) );

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


TEST(GenericType, parameter_during_lazy_finish) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	static Record param_record = Record("record_with_param", "")
			.declare_key("param", Parameter("param"), "desc.")
			.declare_key("some_double", Double(), "Double key")
			.close();

	EXPECT_THROW_WHAT( { TypeBase::lazy_finish(); }, ExcParamaterInIst, "Parameter 'param' appears in the IST");
}


TEST(GenericType, root_without_instance) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	static Record with_root_flag = Record("in_rec", "")
			.root_of_generic_subtree()
			.declare_key("array", Array( Array( Parameter("param") ) ), "desc.")
			.close();

	static Record root_rec = Record("root_record", "desc.")
			.declare_key("primary", with_root_flag, "Primary problem.")
			.declare_key("bool", Bool(), "Some bool key.")
			.close();

	EXPECT_THROW_WHAT( { TypeBase::lazy_finish(); }, ExcGenericWithoutInstance, "'in_rec' used without Instance");
}


template <int spacedim>
class ParametrizedTestClass {
public:
	static string template_name() {
	    stringstream ss;
	    ss << "Rec:" << spacedim;
		return ss.str();
	}

	static Record & get_input_type() {
		static Record rec = Record(template_name(), "Record with parameterized dimension.")
				.root_of_generic_subtree()
				.declare_key("param", Parameter("param"), "desc.")
				.declare_key("boundary", Bool(), "desc.")
				.close();

		return rec;
	}

	static const Instance & get_input_type_instance() {
		std::vector<TypeBase::ParameterPair> param_vec;
		param_vec.push_back( std::make_pair("param", std::make_shared<Array>(Double(), 0, spacedim)) );

		return Instance(get_input_type(), param_vec).close();
	}
};


TEST(GenericType, unused_record) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	{
		// same mechanism as call of registrars
		ParametrizedTestClass<2>::get_input_type();
		ParametrizedTestClass<3>::get_input_type();
	}

	static Record root_templated_record = Record("templated_record", "")
			.declare_key("generic_rec", ParametrizedTestClass<3>::get_input_type_instance(), "desc.")
			.declare_key("some_int", Integer(), "Integer key")
			.close();

	TypeBase::lazy_finish();
}


TEST(GenericType, parameter_not_used) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

	std::vector<TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("param", std::make_shared<Integer>()) );
	param_vec.push_back( std::make_pair("param2", std::make_shared<Double>()) );

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
