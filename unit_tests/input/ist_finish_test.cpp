/**
 * ist_recursion_test.cpp
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include "input/input_type.hh"
#include "input/reader_to_storage.hh"
#include "input/accessors.hh"

namespace IT = Input::Type;

class RecordGeneratorTest {
public:
	static IT::Record &get_root_rec();
	static IT::Abstract &get_abstract();
	static const IT::Record &get_b_rec();
	static const IT::Record &get_c_rec();
	static const IT::Record &get_d_rec();
};

IT::Record &RecordGeneratorTest::get_root_rec() {
    return IT::Record("Root","")
            .declare_key("a_key", RecordGeneratorTest::get_abstract(), "")
            .declare_key("b_key", RecordGeneratorTest::get_b_rec(), "")
			.close();
}

IT::Abstract &RecordGeneratorTest::get_abstract() {
    return IT::Abstract("AbstractWithRecursion","").close();
}

const IT::Record &RecordGeneratorTest::get_b_rec() {
    return IT::Record("B_rec","")
            .declare_key("a_key", RecordGeneratorTest::get_abstract(), "")
            .declare_key("b_val", IT::Integer(),"")
			.close();
}

const IT::Record &RecordGeneratorTest::get_c_rec() {
    return IT::Record("C_rec","")
            .derive_from(RecordGeneratorTest::get_abstract())
            .declare_key("d_key", RecordGeneratorTest::get_d_rec(), "")
            .declare_key("x_val", IT::Integer(),"")
            .declare_key("y_val", IT::Double(),"")
			.close();
}

const IT::Record &RecordGeneratorTest::get_d_rec() {
    return IT::Record("D_rec","")
            .declare_key("a_key", RecordGeneratorTest::get_abstract(), "")
            .declare_key("d_val", IT::Integer(),"")
			.close();
}


TEST(ISTFinish, ist_recursion) {
	EXPECT_EQ( 4, RecordGeneratorTest::get_c_rec().size() ); // touch of record simulates registrar
	IT::Record root = RecordGeneratorTest::get_root_rec();
	EXPECT_ASSERT_DEATH( { root.finish(); }, "AbstractWithRecursion");
}



class DefaultValuesTest {
public:
	static IT::Record &get_root_rec();
	static IT::Abstract &get_abstract();
	static const IT::Record &get_desc();
	static const IT::Record &get_subrecord();
};

IT::Record &DefaultValuesTest::get_root_rec() {
    return IT::Record("RootRecord","")
            .declare_key("data", DefaultValuesTest::get_abstract(), IT::Default("1.0"), "")
            .declare_key("description", IT::String(), "")
			.close();
}

IT::Abstract &DefaultValuesTest::get_abstract() {
    return IT::Abstract("AbstractDefault","")
    	.allow_auto_conversion("DescendantRec")
		.close();
}

const IT::Record &DefaultValuesTest::get_desc() {
    return IT::Record("DescendantRec","")
    		.derive_from(DefaultValuesTest::get_abstract())
            .declare_key("rec_key", DefaultValuesTest::get_subrecord(), IT::Default::obligatory(), "")
            .declare_key("int_val", IT::Integer(), IT::Default("0"), "")
			.allow_auto_conversion("rec_key")
			.close();
}

const IT::Record &DefaultValuesTest::get_subrecord() {
    return IT::Record("SubRecord","")
            .declare_key("val1", IT::Double(), IT::Default::obligatory(), "")
            .declare_key("val2", IT::Integer(), IT::Default("0"), "")
            .declare_key("val3", IT::String(), IT::Default::optional(), "")
			.allow_auto_conversion("val1")
			.close();
}

TEST(ISTFinish, complex_default_val) {
	EXPECT_EQ( 3, DefaultValuesTest::get_desc().size() ); // touch of record simulates registrar
	IT::Record root = DefaultValuesTest::get_root_rec();
	root.finish();
}


class AbstractAddChildTest {
public:
	static IT::Record &get_root_rec();
	static IT::Abstract &get_abstract();
	static const IT::Record &get_desc();
};

IT::Record &AbstractAddChildTest::get_root_rec() {
    return IT::Record("RootRecord","")
            .declare_key("problem", AbstractAddChildTest::get_abstract(), IT::Default::obligatory(), "")
            .declare_key("int_value", IT::Integer(), "")
			.close();
}

IT::Abstract &AbstractAddChildTest::get_abstract() {
	static IT::Abstract abstr = IT::Abstract("Abstract2","")
		.close();
	abstr.add_child( AbstractAddChildTest::get_desc() );

	return abstr;
}

const IT::Record &AbstractAddChildTest::get_desc() {
    return IT::Record("Descendant","")
            .declare_key("int_val", IT::Integer(), IT::Default("0"), "")
            .declare_key("desc", IT::String(), IT::Default::obligatory(), "Description.")
			.close();
}

const string read_input_yaml = R"YAML(
problem: !Descendant
  int_val: 1
  desc: Test case
int_value: 10
)YAML";

TEST(ISTFinish, abstract_add_child) {
	EXPECT_EQ( 3, AbstractAddChildTest::get_desc().size() );

	Input::ReaderToStorage json_reader( read_input_yaml, AbstractAddChildTest::get_root_rec(), Input::FileFormat::format_YAML);
	Input::Record rec=json_reader.get_root_interface<Input::Record>();
	EXPECT_EQ(10, *(rec.find<int>("int_value")) );
	Input::AbstractRecord problem_arec = rec.val<Input::AbstractRecord>("problem");
	Input::Record problem_rec = Input::Record(problem_arec);
	EXPECT_EQ(1, *(Input::Record(problem_arec).find<int>("int_val")) );
	EXPECT_EQ("Test case", *(problem_rec.find<std::string>("desc")) );
	EXPECT_EQ("Descendant", *(problem_rec.find<std::string>("TYPE")) );
}

