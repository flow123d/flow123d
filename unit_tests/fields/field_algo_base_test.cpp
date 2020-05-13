/*
 * field_algo_base_test.cpp
 *
 *  Created on: Feb 3, 2013
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS

#include <memory>
#include <regex>

#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"


#include "fields/field.hh"
#include "fields/field_algo_base.hh"

#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "fields/field_constant.hh"
#include "fields/field_set.hh"
#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"










FLOW123D_FORCE_LINK_IN_PARENT(field_formula)


template <class F>
class FieldFix : public testing::Test, public F {
public:
	typedef F FieldType;
	static constexpr bool is_enum_valued = boost::is_same<typename FieldType::ValueType::element_type, FieldEnum>::value;

	void SetUp() {
	    Profiler::instance();
	    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	    my_mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");


	    field_.name("test_field");
	    field_.input_selection( get_test_selection() );

		auto a_rec_type = this->field_.get_input_type();
		test_field_descriptor = make_shared<Input::Type::Record>(
				Input::Type::Record("any", FieldCommon::field_descriptor_record_description("any") )
				.copy_keys( this->field_.field_descriptor_record("any") )
				.declare_key("a", a_rec_type, "")
				.declare_key("b", a_rec_type, "")
				.close()
				);
		test_input_list = make_shared<Input::Type::Array>( *test_field_descriptor );
		test_input_list->finish();

		my_domain = my_mesh->region_db().get_region_set("BULK");

		string field_input = "[{a=314}]";
		if (this->is_enum_valued) field_input = "[{a=\"white\"}]";

		root_input = input_list(field_input);
		Input::Record x_rec = *(root_input.begin<Input::Record>());
		auto field_rec = *(x_rec.find<Input::AbstractRecord>("a"));
	    FieldAlgoBaseInitData init_data_conc("a", this->n_comp(), UnitSI::dimensionless());
		my_field_algo_base = FieldType::FieldBaseType::function_factory(field_rec, init_data_conc);
	}

	void TearDown() {
		delete my_mesh;
	}

	Input::Array input_list(const string& str) {
		static std::vector<Input::Array> inputs;
		unsigned int input_last = inputs.size(); // position of new item
		Input::ReaderToStorage reader(str, *test_input_list, Input::FileFormat::format_JSON );
		inputs.push_back( reader.get_root_interface<Input::Array>() );
		return inputs[input_last];
	}

	// RegionHistory for given region index
	typename FieldType::RegionHistory &rh( int r_idx) { return this->data_->region_history_[r_idx]; };
	// time of HistoryPoint given by region index and position in circular buffer.
	double rh_time(int r_idx, int j) { return this->rh(r_idx)[j].first; }

	typedef typename FieldType::ValueType Value;
	// const value of HistoryPoint given by region index and position in circular buffer.
	typename Value::element_type rh_value(int r_idx, int j) {
		typename FieldType::FieldBasePtr fb = rh(r_idx)[j].second;
		auto elm = ElementAccessor<FieldType::space_dim>( this->mesh(), 0 );
		auto val = fb->value(elm.centre(), elm);
		return (Value(val))(0,0);
	}
	// const value of region_field_ on some region
	typename Value::element_type _value_(FieldType &f) {
		Region r = f.mesh()->region_db().get_region_set("BULK")[1];

		auto elm = ElementAccessor<FieldType::space_dim>( f.mesh(), r );
		auto val = f.value( point_, elm);
		return (Value(val))(0,0);
	}

	// simple selection with values "black" and "White"
	static const Input::Type::Selection &get_test_selection();

	// has to keep root accessor to prevent delete of the storage tree
	Input::Array root_input;

	// Field<..> instance with name "test_field" and selection test_selection
	FieldType field_;

	/* Regions in the test mesh:
	 * $PhysicalNames
	    6
	    1       37      "1D diagonal"
	    2       38      "2D XY diagonal"
	    2       101     ".top side"
	    2       102     ".bottom side"
	    3       39      "3D back"
	    3       40      "3D front"
	    $EndPhysicalNames
	 */
	Mesh * my_mesh;

	typename FieldType::Point point_;

	// BULK domain
	RegionSet my_domain;

	// FieldConstant with value 314 (numeric return value) or "white" (enum return value)
	std::shared_ptr<typename FieldType::FieldBaseType> my_field_algo_base;

	// simple field descriptor input type conataining fields "a" and "b"
	std::shared_ptr<Input::Type::Record> test_field_descriptor;
	// list of simple field descriptors
	std::shared_ptr<Input::Type::Array>          test_input_list;

};

template <class F>
const Input::Type::Selection &FieldFix<F>::get_test_selection() {
	return Input::Type::Selection("any")
				.add_value(0,"black")
				.add_value(1,"white")
				.close();
}



#define FV FieldValue
// full list
#define f_list(Dim) \
	Field<Dim,FV<0>::Scalar> , \
    Field<Dim,FV<0>::Enum>, \
    Field<Dim,FV<0>::Integer>, \
	Field<Dim,FV<Dim>::VectorFixed>, \
	Field<Dim,FV<Dim>::TensorFixed>

// simple list
#define s_list(Dim) Field<Dim,FV<0>::Scalar>

/**
 * Select set of Value template parameters.
 *
 * Can only use for spacedim==3
 */
typedef ::testing::Types<f_list(3)> FieldTypes;
TYPED_TEST_CASE(FieldFix, FieldTypes);


// check that get_input_type works for every Field<> type
// It should be tested in static context, which is non-trivial,
// we shall do it as part of FieldList test.
//
TYPED_TEST(FieldFix, get_input_type) {
	auto a_rec_type = this->field_.get_input_type();
}


// test correctness of check for ascending time sequence
TYPED_TEST(FieldFix, set_input_list) {
	string list_ok = "["
			"{time=2, a=0}, "
			"{time=1, b=0}, "
			"{time=1, c=0},"
			"{time=3, a=0}, "
			"{time=2, b=0}, "
			"{time=10, c=0},"
			"{time=4, a=0}, "
			"{time=5, a=0, b=0}]";

	string list_ko = "["
			"{time=2, a=0},"
			"{time=1, a=0},"
			"{time=2, a=0}]";

	if (this->is_enum_valued) {
		std::regex e("=0");
		list_ok = std::regex_replace(list_ok, e, "=\"white\"");
		list_ko = std::regex_replace(list_ko, e, "=\"white\"");
	}

	TimeGovernor tg(0.0, 1.0);
	this->field_.name("a");
	this->field_.set_input_list( this->input_list(list_ok), tg );

	this->field_.name("b");
	this->field_.set_input_list( this->input_list(list_ok), tg );

	this->field_.name("a");
	EXPECT_THROW_WHAT( {	this->field_.set_input_list( this->input_list(list_ko), tg );}, FieldCommon::ExcNonascendingTime, "for field 'a'" );
}


// check that it correctly introduce requered marks
TYPED_TEST(FieldFix, mark_input_times) {
    string list_ok = "["
            "{time=2, a=0}, "
            "{time=1, b=0}, "
            "{time=1, c=0},"
            "{time=3, a=0}, "
            "{time=2, b=0}, "
            "{time=10, c=0},"
            "{time=4, a=0}, "
            "{time=5, a=0, b=0}]";

	if (this->is_enum_valued) {
		std::regex e("=0");
		list_ok = std::regex_replace(list_ok, e, "=\"white\"");
	}

	TimeGovernor tg;

	this->field_.name("b");
	this->field_.set_input_list(this->input_list(list_ok), tg);

	TimeMark::Type mark_type = TimeMark::Type(tg.marks().type_fixed_time().bitmap_, tg.equation_mark_type().equation_index_);
	this->field_.mark_input_times(tg);
	auto it = tg.marks().next(tg, mark_type);
	EXPECT_EQ( 1, it->time());
	EXPECT_TRUE( it->match_mask(mark_type) );
	++it;
	EXPECT_EQ( 2, it->time());
	EXPECT_TRUE( it->match_mask(mark_type) );
	++it;
	EXPECT_EQ( 5, it->time());
	EXPECT_TRUE( it->match_mask(mark_type) );
	++it;
	EXPECT_EQ( TimeGovernor::inf_time, it->time());
	EXPECT_TRUE( it->match_mask(mark_type) );

}


TYPED_TEST(FieldFix, set_mesh) {
	EXPECT_ASSERT_DEATH( {this->set_field(this->my_domain, this->my_field_algo_base);}, "Null; mesh pointer");

	EXPECT_EQ(nullptr, this->shared_->mesh_);
	this->set_mesh(*(this->my_mesh));
	EXPECT_EQ(this->my_mesh, this->shared_->mesh_);
	int size = this->my_mesh->region_db().size();
	EXPECT_EQ( size, this->region_fields_.size());
	EXPECT_EQ( size, this->data_->region_history_.size());


}



TYPED_TEST(FieldFix, set_field) {
	this->set_mesh(*(this->my_mesh));
	this->set_field(this->my_domain, this->my_field_algo_base);

	Region reg = this->my_domain[0];
	auto const &history = this->data_->region_history_[reg.idx()];

	EXPECT_EQ(1, history.size());
	EXPECT_EQ(0.0, history[0].first);
	EXPECT_TRUE( bool(history[0].second) );

	this->set_field(this->my_domain, this->my_field_algo_base, 3.0);
	EXPECT_EQ(2, history.size());
	EXPECT_EQ(3.0, history[0].first);
	EXPECT_EQ(0.0, history[1].first);

	EXPECT_ASSERT_DEATH( {this->set_field(this->my_domain, this->my_field_algo_base, 1.0);}, "" );

	this->set_field(this->my_domain, this->my_field_algo_base, 6.0);
	EXPECT_EQ(3, history.size());
	EXPECT_EQ(6.0, history[0].first);
	EXPECT_EQ(3.0, history[1].first);
	EXPECT_EQ(0.0, history[2].first);

	this->set_field(this->my_domain, this->my_field_algo_base, 7.0);
	EXPECT_EQ(3, history.size());
	EXPECT_EQ(7.0, history[0].first);
	EXPECT_EQ(6.0, history[1].first);
	EXPECT_EQ(3.0, history[2].first);

}


TYPED_TEST(FieldFix, update_history) {
	string list_ok = "["
			"{time=0, region=\"ALL\", a =0, b =0},"
			"{time=1, region=\"BULK\", a =1, b =0},"
			"{time=2, region=\".BOUNDARY\", a =1, b =0},"
			"{time=3, region=\"ALL\", b =0},"
			"{time=4, region=\"ALL\", a =0},"
			"{time=5, region=\"ALL\", a =1}"
			"]";
	if (this->is_enum_valued) {
		list_ok = std::regex_replace(list_ok, std::regex(" =1"), "=\"white\"");
		list_ok = std::regex_replace(list_ok, std::regex(" =0"), "=\"black\"");
	}

	TimeGovernor tg(0.0, 1.0);
	this->name("a");
	this->set_mesh(*(this->my_mesh));
	this->set_input_list( this->input_list(list_ok), tg );
	this->units( UnitSI().m() );

	// time = 0.0
	this->update_history(tg.step());

    Region diagonal_1d = this->mesh()->region_db().find_label("1D diagonal");
    Region diagonal_2d = this->mesh()->region_db().find_label("2D XY diagonal");
    Region bc_top = this->mesh()->region_db().find_label(".top side");
    Region bc_bottom = this->mesh()->region_db().find_label(".bottom side");
    Region front_3d = this->mesh()->region_db().find_label("3D front");
    Region back_3d = this->mesh()->region_db().find_label("3D back");

    EXPECT_EQ( 1 , this->rh(diagonal_1d.idx()).size() );
	EXPECT_EQ( 1 , this->rh(diagonal_2d.idx()).size() );
	EXPECT_EQ( 1 , this->rh(bc_top.idx()).size() );
	EXPECT_EQ( 1 , this->rh(bc_bottom.idx()).size() );
	EXPECT_EQ( 1 , this->rh(front_3d.idx()).size() );
	EXPECT_EQ( 1 , this->rh(back_3d.idx()).size() );

	EXPECT_EQ( 0.0 , this->rh_time(front_3d.idx(),0) );
	EXPECT_EQ( 0.0 , this->rh_time(bc_top.idx(),0) );

	EXPECT_EQ( 0 , this->rh_value(front_3d.idx(),0) );
	EXPECT_EQ( 0 , this->rh_value(bc_top.idx(),0) );

	tg.estimate_dt();
	tg.next_time();
	this->update_history(tg.step());

	// time = 1.0
    EXPECT_EQ( 2 , this->rh(diagonal_1d.idx()).size() );
	EXPECT_EQ( 2 , this->rh(diagonal_2d.idx()).size() );
	EXPECT_EQ( 1 , this->rh(bc_top.idx()).size() );
	EXPECT_EQ( 1 , this->rh(bc_bottom.idx()).size() );
	EXPECT_EQ( 2 , this->rh(front_3d.idx()).size() );
	EXPECT_EQ( 2 , this->rh(back_3d.idx()).size() );

	EXPECT_EQ( 1.0 , this->rh_time(front_3d.idx(),0) );
	EXPECT_EQ( 0.0 , this->rh_time(bc_top.idx(),0) );

	EXPECT_EQ( 1 , this->rh_value(front_3d.idx(),0) );
	EXPECT_EQ( 0 , this->rh_value(bc_top.idx(),0) );

	tg.next_time();
	this->update_history(tg.step());

	// time = 2.0
    EXPECT_EQ( 2 , this->rh(diagonal_1d.idx()).size() );
	EXPECT_EQ( 2 , this->rh(diagonal_2d.idx()).size() );
	EXPECT_EQ( 2 , this->rh(bc_top.idx()).size() );
	EXPECT_EQ( 2 , this->rh(bc_bottom.idx()).size() );
	EXPECT_EQ( 2 , this->rh(front_3d.idx()).size() );
	EXPECT_EQ( 2 , this->rh(back_3d.idx()).size() );

	EXPECT_EQ( 1.0 , this->rh_time(front_3d.idx(),0) );
	EXPECT_EQ( 2.0 , this->rh_time(bc_top.idx(),0) );

	EXPECT_EQ( 1 , this->rh_value(front_3d.idx(),0) );
	EXPECT_EQ( 1 , this->rh_value(bc_top.idx(),0) );

	tg.next_time();
	this->update_history(tg.step());

	// time = 3.0
    EXPECT_EQ( 2 , this->rh(diagonal_1d.idx()).size() );
	EXPECT_EQ( 2 , this->rh(diagonal_2d.idx()).size() );
	EXPECT_EQ( 2 , this->rh(bc_top.idx()).size() );
	EXPECT_EQ( 2 , this->rh(bc_bottom.idx()).size() );
	EXPECT_EQ( 2 , this->rh(front_3d.idx()).size() );
	EXPECT_EQ( 2 , this->rh(back_3d.idx()).size() );

	EXPECT_EQ( 1.0 , this->rh_time(front_3d.idx(),0) );
	EXPECT_EQ( 2.0 , this->rh_time(bc_top.idx(),0) );

	EXPECT_EQ( 1 , this->rh_value(front_3d.idx(),0) );
	EXPECT_EQ( 1 , this->rh_value(bc_top.idx(),0) );

	tg.next_time();
	this->update_history(tg.step());

	// time = 4.0
    EXPECT_EQ( 3 , this->rh(diagonal_1d.idx()).size() );
	EXPECT_EQ( 3 , this->rh(diagonal_2d.idx()).size() );
	EXPECT_EQ( 3 , this->rh(bc_top.idx()).size() );
	EXPECT_EQ( 3 , this->rh(bc_bottom.idx()).size() );
	EXPECT_EQ( 3 , this->rh(front_3d.idx()).size() );
	EXPECT_EQ( 3 , this->rh(back_3d.idx()).size() );

	EXPECT_EQ( 4.0 , this->rh_time(front_3d.idx(),0) );
	EXPECT_EQ( 4.0 , this->rh_time(bc_top.idx(),0) );

	EXPECT_EQ( 0 , this->rh_value(front_3d.idx(),0) );
	EXPECT_EQ( 0 , this->rh_value(bc_top.idx(),0) );

	tg.next_time();
	this->update_history(tg.step());

	// time = 5.0
    EXPECT_EQ( 3 , this->rh(diagonal_1d.idx()).size() );
	EXPECT_EQ( 3 , this->rh(diagonal_2d.idx()).size() );
	EXPECT_EQ( 3 , this->rh(bc_top.idx()).size() );
	EXPECT_EQ( 3 , this->rh(bc_bottom.idx()).size() );
	EXPECT_EQ( 3 , this->rh(front_3d.idx()).size() );
	EXPECT_EQ( 3 , this->rh(back_3d.idx()).size() );

	EXPECT_EQ( 5.0 , this->rh_time(front_3d.idx(),0) );
	EXPECT_EQ( 5.0 , this->rh_time(bc_top.idx(),0) );

	EXPECT_EQ( 1 , this->rh_value(front_3d.idx(),0) );
	EXPECT_EQ( 1 , this->rh_value(bc_top.idx(),0) );

}

TYPED_TEST(FieldFix, set_time) {
	string list_ok = "["
			"{time=0, region=\"ALL\", a =0, b =0},"
			"{time=1, region=\"BULK\", a =1, b =0},"
			"{time=2, region=\".BOUNDARY\", a =1, b =0},"
			"{time=3, region=\"ALL\", b =0},"
			"{time=4, region=\"ALL\", a =0},"
			"{time=5, region=\"ALL\", a =1}"
			"]";

	if (this->is_enum_valued) {
		list_ok = std::regex_replace(list_ok, std::regex(" =1"), "=\"white\"");
		list_ok = std::regex_replace(list_ok, std::regex(" =0"), "=\"black\"");
	}

	TimeGovernor tg(0.0, 0.5);
	this->name("a");
	this->set_mesh(*(this->my_mesh));
	this->set_input_list( this->input_list(list_ok), tg );
	this->units( UnitSI().m() );

	// time = 0.0
	this->set_time(tg.step(), LimitSide::right);
	EXPECT_EQ(0, this->_value_( *this ));
	EXPECT_TRUE( this->is_jump_time() );

	tg.next_time();
	this->set_time(tg.step(), LimitSide::left);
	EXPECT_EQ(0, this->_value_( *this ));
    EXPECT_FALSE( this->is_jump_time() );

    this->set_time(tg.step(), LimitSide::right);
    EXPECT_EQ(0, this->_value_( *this ));
    EXPECT_FALSE( this->is_jump_time() );

    tg.next_time();
    this->set_time(tg.step(), LimitSide::left);
    EXPECT_EQ(0, this->_value_( *this ));
    EXPECT_TRUE( this->is_jump_time() );

    this->set_time(tg.step(), LimitSide::right);
    EXPECT_EQ(1, this->_value_( *this ));
    EXPECT_TRUE( this->is_jump_time() );

}



// Check copy constructor and assignment oprerator
TYPED_TEST(FieldFix, constructors) {
	// default constructor
	typename TestFixture::FieldType field_default;
	EXPECT_EQ("", field_default.name());
	EXPECT_FALSE(field_default.is_bc());

	// copies
	// check that we can have copies in different times
	this->field_.name("a");
	field_default
	    .name("b")
	    .flags(FieldFlag::input_copy);
	this->field_.set_mesh( *(this->my_mesh) );
	field_default.set_mesh( *(this->my_mesh) );
	this->field_.units( UnitSI().m() );
	field_default.units( UnitSI().m() );

	string list_ok = "["
			"{time=2,  region=\"BULK\", a=0, b=1}, "
			"{time=3,  region=\"BULK\", b=1}, "
			"{time=4,  region=\"BULK\", a=1},"
			"{time=5,  region=\"BULK\", a=0, b=0}]";

	if (this->is_enum_valued) {
		list_ok = std::regex_replace(list_ok, std::regex("=1"), "=\"white\"");
		list_ok = std::regex_replace(list_ok, std::regex("=0"), "=\"black\"");
	}

	TimeGovernor tg(2.0, 1.0);
	this->field_.set_input_list(this->input_list(list_ok), tg);
	field_default.set_input_list(this->input_list(list_ok), tg);



	typename TestFixture::FieldType f2(this->field_);	// default constructor
	field_default = this->field_; // assignment, should overwrite name "b" by name "a"


	// tg = 2.0
	f2.set_time(tg.step(), LimitSide::right);
	EXPECT_EQ(0,this->_value_(f2));
	EXPECT_ASSERT_DEATH( {this->_value_(this->field_);}, "");
	this->field_.set_time(tg.step(), LimitSide::right);
	EXPECT_EQ(0,this->_value_(this->field_));

	// tg = 3.0
	tg.next_time();
	this->field_.set_time(tg.step(), LimitSide::right);
	EXPECT_EQ(0,this->_value_(this->field_));

	// tg = 4.0
	tg.next_time();
	this->field_.set_time(tg.step(), LimitSide::right);
	EXPECT_EQ(1,this->_value_(this->field_));
	EXPECT_EQ(0,this->_value_(f2));

	field_default.set_time(tg.step(), LimitSide::right);
	EXPECT_EQ(1,this->_value_(field_default));
}







string field_input = R"INPUT(
{
   sorption_type="linear",   
   init_conc=[ 10, 20, 30],    // FieldConst
   conductivity={ //3x3 tensor
       TYPE="FieldFormula",
       value=["x","y", "z"]
   },
   conductivity_3d={ //3x3 tensor - for test of Field::is_constant method
       TYPE="FieldFormula",
       value=["1", "t", "t*t"]
   }
}
)INPUT";

namespace it = Input::Type;

static const it::Selection &get_sorption_type_selection() {
	return it::Selection("SorptionType")
				.add_value(1,"linear")
				.add_value(0,"none")
				.close();
}

TEST(Field, init_from_input) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";
	Profiler::instance();
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");

    Field<3, FieldValue<3>::Enum > sorption_type;
    Field<3, FieldValue<3>::VectorFixed > init_conc;
    Field<3, FieldValue<3>::TensorFixed > conductivity;
    Field<3, FieldValue<3>::TensorFixed > conductivity_3d;


    std::vector<string> component_names = { "comp_0", "comp_1", "comp_2" };


    sorption_type.input_selection( get_sorption_type_selection() );
    init_conc.set_components(component_names);

    it::Record main_record =
            it::Record("main", "desc")
            .declare_key("sorption_type", sorption_type.get_input_type(), it::Default::obligatory(), "desc")
            .declare_key("init_conc", init_conc.get_input_type(), it::Default::obligatory(), "desc")
            .declare_key("conductivity", conductivity.get_input_type(), it::Default::obligatory(), "desc")
            .declare_key("conductivity_3d", conductivity_3d.get_input_type(), it::Default::obligatory(), "desc")
			.close();


    // read input string
    Input::ReaderToStorage reader( field_input, main_record, Input::FileFormat::format_JSON );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    sorption_type.set_mesh(*mesh);
    init_conc.set_mesh(*mesh);
    conductivity.set_mesh(*mesh);
    conductivity_3d.set_mesh(*mesh);

    sorption_type.units( UnitSI().m() );
    init_conc.units( UnitSI().m() );
    conductivity.units( UnitSI().m() );
    conductivity_3d.units( UnitSI().m() );

    auto region_set = mesh->region_db().get_region_set("BULK");

    sorption_type.set_field(region_set, in_rec.val<Input::AbstractRecord>("sorption_type"));
    init_conc.set_field(region_set, in_rec.val<Input::AbstractRecord>("init_conc"));
    conductivity.set_field(region_set, in_rec.val<Input::AbstractRecord>("conductivity"));
    conductivity_3d.set_field(region_set, in_rec.val<Input::AbstractRecord>("conductivity_3d"));





    sorption_type.set_time(TimeGovernor().step(), LimitSide::right);
    init_conc.set_time(TimeGovernor().step(), LimitSide::right);
    conductivity.set_time(TimeGovernor().step(), LimitSide::right);
    conductivity_3d.set_time(TimeGovernor().step(), LimitSide::right);

    {	

	    auto ele = mesh->element_accessor(5);

	    EXPECT_EQ( 1, sorption_type.value(ele.centre(), ele) );


	    auto vec_value = init_conc.value(ele.centre(), ele);
	    EXPECT_TRUE( arma::min( arma::vec("10 20 30") == vec_value ) );

	    auto result =conductivity.value(ele.centre(), ele);
	    arma::mat diff = arma::mat33("-0.5 0 0;0 0 0; 0 0 -0.5") - result;

	    double norm=arma::norm(diff, 1);
	    EXPECT_DOUBLE_EQ( 0.0, norm );
    }

    {
	//  using const accessor
    	ElementAccessor<3> ele;

	    EXPECT_ASSERT_DEATH( {sorption_type.value(ele.centre(), ele);}  , "Invalid element accessor.");
	    Region reg = mesh->region_db().find_id(40);

	    EXPECT_TRUE( sorption_type.is_constant(reg) );
	    EXPECT_TRUE( init_conc.is_constant(reg) );
	    EXPECT_FALSE( conductivity.is_constant(reg) );
	    EXPECT_TRUE( conductivity_3d.is_constant(reg) );

	    ele = ElementAccessor<3>(mesh, reg);
	    EXPECT_EQ( 1, sorption_type.value(ele.centre(), ele) );
	    EXPECT_TRUE( arma::min( arma::vec("10 20 30") == init_conc.value(ele.centre(), ele) ) );

   }

    delete mesh;
}



string field_input_list = R"INPUT(
[
    {
        region="1D diagonal",
        scalar=0,
        vector=0,
        tensor=0
    },
    {
        region="2D XY diagonal",
        scalar=1,
        vector=1,
        tensor=[1,1,1,1,1,1]
    },
    {
        region="3D front",
        scalar=2,
        vector=2,
        tensor=2
    },
    {
        region="3D back",
        scalar={TYPE="FieldFormula", value="0"},
        vector={TYPE="FieldFormula", value="0"},
        tensor=1
    }
]
)INPUT";

namespace it = Input::Type;

class TestFieldSet : public FieldSet
{
public:
    TestFieldSet() {
        *this += scalar.name("scalar").units(UnitSI::dimensionless());
        *this += vector.name("vector").units(UnitSI::dimensionless());
        *this += tensor.name("tensor").units(UnitSI::dimensionless());
    }
    Field<3, FieldValue<3>::Scalar > scalar;
    Field<3, FieldValue<3>::VectorFixed > vector;
    Field<3, FieldValue<3>::TensorFixed > tensor;
};

TEST(Field, field_result) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Profiler::instance();
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    TimeGovernor tg(0.0, 1.0);

    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");

    it::Array main_array =IT::Array(
            TestFieldSet().make_field_descriptor_type("TestFieldSet")
            .close()
        );

    // read input string
    Input::ReaderToStorage reader( field_input_list, main_array, Input::FileFormat::format_JSON );
    Input::Array array=reader.get_root_interface<Input::Array>();

    TestFieldSet data;
    data.set_mesh(*mesh);
    data.set_input_list(array, tg);


    Region diagonal_1d = mesh->region_db().find_label("1D diagonal");
    Region diagonal_2d = mesh->region_db().find_label("2D XY diagonal");
    Region front_3d = mesh->region_db().find_label("3D front");
    Region back_3d = mesh->region_db().find_label("3D back");
    Region top_side = mesh->region_db().find_label(".top side");
    Region bottom_side = mesh->region_db().find_label(".bottom side");

    EXPECT_EQ( result_none, data.scalar.field_result({diagonal_1d}) );
    EXPECT_EQ( result_none, data.scalar.field_result({diagonal_2d}) );
    EXPECT_EQ( result_none, data.vector.field_result({diagonal_1d}) );
    EXPECT_EQ( result_none, data.vector.field_result({front_3d}) );
    EXPECT_EQ( result_none, data.tensor.field_result({diagonal_1d}) );
    EXPECT_EQ( result_none, data.tensor.field_result({back_3d}) );
    // time 0
    data.set_time(tg.step(), LimitSide::right);
    EXPECT_EQ( result_zeros, data.scalar.field_result({diagonal_1d}) );
    EXPECT_EQ( result_zeros, data.vector.field_result({diagonal_1d}) );
    EXPECT_EQ( result_zeros, data.tensor.field_result({diagonal_1d}) );

    EXPECT_EQ( result_ones, data.scalar.field_result({diagonal_2d}) );
    EXPECT_EQ( result_ones, data.vector.field_result({diagonal_2d}) );
    EXPECT_EQ( result_ones, data.tensor.field_result({diagonal_2d}) );

    EXPECT_EQ( result_constant, data.tensor.field_result({front_3d}) );
    EXPECT_EQ( result_constant, data.tensor.field_result({front_3d}) );
    EXPECT_EQ( result_constant, data.tensor.field_result({front_3d}) );

    EXPECT_EQ( result_other, data.scalar.field_result({back_3d}) );
    EXPECT_EQ( result_other, data.vector.field_result({back_3d}) );
    EXPECT_EQ( result_eye, data.tensor.field_result({back_3d}) );


    EXPECT_EQ( result_other, data.scalar.field_result({diagonal_1d, diagonal_2d}) );
    EXPECT_EQ( result_other, data.vector.field_result({diagonal_1d, diagonal_2d}) );
    EXPECT_EQ( result_other, data.tensor.field_result({diagonal_1d, diagonal_2d}) );

    delete mesh;
}







static const it::Selection &get_test_type_selection() {
	return it::Selection("TestType")
				.add_value(0, "none")
				.add_value(1,"dirichlet")
				.close();
}

TEST(Field, init_from_default) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Profiler::instance();
    
    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");

    Space<3>::Point p("1 2 3");

    {
        Field<3, FieldValue<3>::Scalar > scalar_field("scalar_test");

        // test default initialization of scalar field
        scalar_field.input_default( "45.0" );
        scalar_field.set_mesh(*mesh);
        scalar_field.units( UnitSI().m() );

        scalar_field.set_time(TimeGovernor().step(), LimitSide::right);

        EXPECT_EQ( 45.0, scalar_field.value(p, mesh->element_accessor(0)) );
        EXPECT_EQ( 45.0, scalar_field.value(p, mesh->element_accessor(6)) );
        // this fails on dev.nti.tul.cz
        //EXPECT_DEATH( { scalar_field.value(p, mesh->element_accessor(0,true)); }, "Null field ptr " );
    }

    {
        Field<3, FieldValue<3>::Scalar > scalar_field("some", true);

        // test death of set_time without default value
        scalar_field.set_mesh(*mesh);

        EXPECT_THROW_WHAT( {scalar_field.set_time(TimeGovernor().step(), LimitSide::right);} , ExcXprintfMsg, "Missing value of the input field");
    }

    {
        Field<3, FieldValue<3>::Enum > enum_field("any", true);

        enum_field.input_selection( get_test_type_selection() );
        enum_field.input_default( "\"none\"" );
        enum_field.set_mesh(*mesh);
        enum_field.units( UnitSI().m() );

        enum_field.set_time(TimeGovernor().step(), LimitSide::right);

        EXPECT_EQ( 0 , enum_field.value(p, mesh->element_accessor(12)) );

    }

    delete mesh;
}

/// Test optional fields dependent e.g. on BC type
TEST(Field, disable_where) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

    enum {
        dirichlet,
        neumann,
        robin
    };
    // test optional checking in the set_time method
    Field<3, FieldValue<3>::Enum > bc_type("bc_type", true);

    std::vector<FieldEnum> list;
    Field<3, FieldValue<3>::Scalar > bc_value("bc_value", true);
    bc_value.disable_where( bc_type, {neumann} );

    Field<3, FieldValue<3>::Scalar > bc_flux("bc_flux", true);
    bc_flux.disable_where( bc_type, {dirichlet, robin} );

    Field<3, FieldValue<3>::Scalar > bc_sigma("bc_sigma", true);
    bc_sigma.disable_where( bc_type, {dirichlet, neumann} );

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");

    bc_type.set_mesh(*mesh);
    bc_flux.set_mesh(*mesh);
    bc_value.set_mesh(*mesh);
    bc_sigma.set_mesh(*mesh);

    /*
    1       37      "1D diagonal"
    2       38      "2D XY diagonal"
    2       101     ".top side"
    2       102     ".bottom side"
    3       39      "3D back"
    3       40      "3D front"
     */


    typedef FieldConstant<3, FieldValue<3>::Scalar > SConst;
    typedef FieldConstant<3, FieldValue<3>::Enum > EConst;
    auto neumann_type = std::make_shared<EConst>();
    neumann_type->set_value(neumann);
    auto robin_type = std::make_shared<EConst>();
    robin_type->set_value(robin);
    auto one = std::make_shared<SConst>();
    one->set_value(1.0);

    bc_type.set_field(RegionSet(1, mesh->region_db().find_id(101)), neumann_type );
    bc_flux.set_field(RegionSet(1, mesh->region_db().find_id(101)), one );

    bc_type.set_field(RegionSet(1, mesh->region_db().find_id(102)), robin_type );
    bc_value.set_field(RegionSet(1, mesh->region_db().find_id(102)), one );
    bc_sigma.set_field(RegionSet(1, mesh->region_db().find_id(102)), one );

    bc_type.set_field(RegionSet(1, mesh->region_db().find_id(-3)), neumann_type );
    bc_flux.set_field(RegionSet(1, mesh->region_db().find_id(-3)), one );






    bc_type.set_time(TimeGovernor().step(), LimitSide::right);
    bc_flux.set_time(TimeGovernor().step(), LimitSide::right);
    bc_value.set_time(TimeGovernor().step(), LimitSide::right);
    bc_sigma.set_time(TimeGovernor().step(), LimitSide::right);

    delete mesh;
}



string field_value_input = R"INPUT(
{
    color="blue",   
    integer=-1,
    scalar=1.5,
    vector=[1, 2, 3],
    tensor=[4, 5, 6]
}
)INPUT";


static const it::Selection &get_color_selection() {
	return it::Selection("ColorType")
				.add_value(1,"blue")
				.add_value(0,"red")
				.close();
}

TEST(Field, field_values) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";
	Profiler::instance();
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/cube_2x1.msh\"}");
	std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);

    Field<3, FieldValue<0>::Enum > color_field;
    Field<3, FieldValue<0>::Integer > int_field;
    Field<3, FieldValue<3>::Scalar > scalar_field;
    Field<3, FieldValue<3>::VectorFixed > vector_field;
    Field<3, FieldValue<3>::TensorFixed > tensor_field;


    //std::vector<string> component_names = { "comp_0", "comp_1", "comp_2" };
    color_field.input_selection( get_color_selection() );
    //init_conc.set_components(component_names);

    it::Record main_record =
            it::Record("main", "desc")
            .declare_key("color", color_field.get_input_type(), it::Default::obligatory(), "desc")
            .declare_key("integer", int_field.get_input_type(), it::Default::obligatory(), "desc")
            .declare_key("scalar", scalar_field.get_input_type(), it::Default::obligatory(), "desc")
            .declare_key("vector", vector_field.get_input_type(), it::Default::obligatory(), "desc")
            .declare_key("tensor", tensor_field.get_input_type(), it::Default::obligatory(), "desc")
			.close();


    // read input string
    Input::ReaderToStorage reader( field_value_input, main_record, Input::FileFormat::format_JSON );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    color_field.set_mesh(*mesh);
    int_field.set_mesh(*mesh);
    scalar_field.set_mesh(*mesh);
    vector_field.set_mesh(*mesh);
    tensor_field.set_mesh(*mesh);

    color_field.units( UnitSI().m() );
    int_field.units( UnitSI().m() );
    scalar_field.units( UnitSI().m() );
    vector_field.units( UnitSI().m() );
    tensor_field.units( UnitSI().m() );

    auto region_set = mesh->region_db().get_region_set("BULK");

    color_field.set_field(region_set, in_rec.val<Input::AbstractRecord>("color"));
    int_field.set_field(region_set, in_rec.val<Input::AbstractRecord>("integer"));
    scalar_field.set_field(region_set, in_rec.val<Input::AbstractRecord>("scalar"));
    vector_field.set_field(region_set, in_rec.val<Input::AbstractRecord>("vector"));
    tensor_field.set_field(region_set, in_rec.val<Input::AbstractRecord>("tensor"));

    color_field.set_time(TimeGovernor().step(), LimitSide::right);
    int_field.set_time(TimeGovernor().step(), LimitSide::right);
    scalar_field.set_time(TimeGovernor().step(), LimitSide::right);
    vector_field.set_time(TimeGovernor().step(), LimitSide::right);
    tensor_field.set_time(TimeGovernor().step(), LimitSide::right);

    // initialize and allocate FieldValueCaches
    std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    std::shared_ptr<BulkIntegral> mass_eval = eval_points->add_bulk<3>(*q_bulk );
    std::shared_ptr<EdgeIntegral> side_eval = eval_points->add_edge<3>(*q_side );
    ElementCacheMap elm_cache_map;
    elm_cache_map.init(eval_points);
    color_field.cache_allocate(eval_points);
    int_field.cache_allocate(eval_points);
    scalar_field.cache_allocate(eval_points);
    vector_field.cache_allocate(eval_points);
    tensor_field.cache_allocate(eval_points);

    // fill FieldValueCaches
    DHCellAccessor dh_cell(dh.get(), 4);
    elm_cache_map.start_elements_update();
    elm_cache_map.add(dh_cell);
    elm_cache_map.prepare_elements_to_update();
    elm_cache_map.mark_used_eval_points( dh_cell, mass_eval->get_subset_idx(), eval_points->subset_size(dh_cell.dim(), mass_eval->get_subset_idx()) );
    elm_cache_map.create_elements_points_map();
    color_field.cache_update(elm_cache_map);
    int_field.cache_update(elm_cache_map);
    scalar_field.cache_update(elm_cache_map);
    vector_field.cache_update(elm_cache_map);
    tensor_field.cache_update(elm_cache_map);
    elm_cache_map.finish_elements_update();

    DHCellAccessor cache_cell = elm_cache_map(dh_cell);
    for(BulkPoint q_point: mass_eval->points(cache_cell, &elm_cache_map)) {
        EXPECT_EQ( 1, color_field(q_point) );
        EXPECT_EQ( -1, int_field(q_point) );
        EXPECT_DOUBLE_EQ( 1.5, scalar_field(q_point) );
        EXPECT_ARMA_EQ( arma::vec3("1 2 3"), vector_field(q_point) );
        EXPECT_ARMA_EQ( arma::mat33("4 0 0; 0 5 0; 0 0 6"), tensor_field(q_point) );
    }

    delete mesh;
}



