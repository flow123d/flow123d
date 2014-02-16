/*
 * field_base_test.cpp
 *
 *  Created on: Feb 3, 2013
 *      Author: jb
 */

#define TEST_USE_MPI

#include <memory>
#include <boost/regex.hpp>

#include <flow_gtest_mpi.hh>


#include "fields/field.hh"
#include "fields/field_base.hh"

#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"
#include "fields/field_constant.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"



string field_input = R"INPUT(
{
   sorption_type="linear",   
   init_conc=[ 10, 20, 30],    // FieldConst
   conductivity={ //3x3 tensor
       TYPE="FieldFormula",
       value=["x","y", "z"]
   }
}
)INPUT";





template <class F>
class FieldFix : public testing::Test, public F {
public:
	typedef F FieldType;

	void SetUp() {
	    Profiler::initialize();

		test_selection = Input::Type::Selection("any")
			.add_value(0,"black")
			.add_value(1,"white");

	    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh", FilePath::input_file);
	    my_mesh = new Mesh();
	    ifstream in(string(mesh_file).c_str());
	    my_mesh->read_gmsh_from_stream(in);


	    field_.name("test_field");
	    field_.input_selection(test_selection);

		auto a_rec_type = this->field_.get_input_type();
		test_field_descriptor = make_shared<Input::Type::Record>(
				this->field_.field_descriptor_record("any")
				.declare_key("a", a_rec_type, "")
				.declare_key("b", a_rec_type, "")
				);
		test_input_list = make_shared<Input::Type::Array>( *test_field_descriptor );
		test_input_list->finish();

		my_domain = my_mesh->region_db().get_region_set("BULK");

		string field_input = "[{a=314}]";
		if (this->is_enum_valued) field_input = "[{a=\"white\"}]";
		Input::Record x_rec = *(input_list(field_input)
						.begin<Input::Record>());
		auto field_rec = *(x_rec.find<Input::AbstractRecord>("a"));
		my_field_base = FieldType::FieldBaseType::function_factory(field_rec, this->n_comp());

		//input_list("[]");
	}

	void TearDown() {
		delete my_mesh;
	}

	Input::Array input_list(const string& str) {
		reader_.read_from_string( str, *test_input_list );
		return reader_.get_root_interface<Input::Array>();
	}

	typename FieldType::RegionHistory &rh( int r_idx) { return this->data_->region_history_[r_idx]; };
	double rh_time(int r_idx, int j) { return this->rh(r_idx)[j].first; }

	typedef typename FieldType::FieldBaseType::ValueType Value;
	typename Value::element_type rh_value(int r_idx, int j) {
		typename FieldType::FieldBasePtr fb = rh(r_idx)[j].second;
		auto elm = ElementAccessor<FieldType::space_dim>( this->mesh(), 0 , this->is_bc() );
		auto val = fb->value(elm.centre(), elm);
		return (Value(val))(0,0);
	}

	Input::JSONToStorage reader_;

	// simple selection with values "black" and "White"
	Input::Type::Selection test_selection;

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
	Mesh *my_mesh;

	// BULK domain
	RegionSet my_domain;

	// FieldConstant with value 314 (numeric return value) or "white" (enum return value)
	std::shared_ptr<typename FieldType::FieldBaseType> my_field_base;

	// simple field descriptor input type conataining fields "a" and "b"
	std::shared_ptr<Input::Type::Record> test_field_descriptor;
	// list of simple field descriptors
	std::shared_ptr<Input::Type::Array>          test_input_list;

};



#define FV FieldValue
#define f_list(Dim) \
	Field<Dim,FV<0>::Scalar> , \
	Field<Dim,FV<0>::Vector>, \
    Field<Dim,FV<0>::Enum>, \
    Field<Dim,FV<0>::EnumVector>, \
    Field<Dim,FV<0>::Integer>, \
	Field<Dim,FV<0>::Vector>, \
	Field<Dim,FV<2>::VectorFixed>, \
	Field<Dim,FV<3>::VectorFixed>, \
	Field<Dim,FV<2>::TensorFixed>, \
	Field<Dim,FV<3>::TensorFixed>

#define s_list(Dim) Field<Dim,FV<0>::Scalar>

/**
 * Select set of Value template parameters.
 *
 * Can only use for spacedim==3
 */
//typedef ::testing::Types< s_list(3)> FieldTypes;
typedef ::testing::Types<f_list(3)> FieldTypes;
TYPED_TEST_CASE(FieldFix, FieldTypes);


// check that get_input_type works for every Field<> type
// It should be tested in static context, which is non-trivial,
// we shall do it as part of FieldList test.
//
TYPED_TEST(FieldFix, get_input_type) {
	Input::Type::AbstractRecord a_rec_type = this->field_.get_input_type();
}


// test correctness of check for ascending time sequence
TYPED_TEST(FieldFix, set_input_list) {
	string list_ok = "["
			"{time=2, a=0}, {time=1, b=0}, {time=10, c=0},"
			"{time=3, a=0}, {time=2, b=0}, {time=1, c=0},"
			"{time=4, a=0}, {time=5, a=0, b=0}]";
	string list_ko = "["
			"{time=2, a=0},"
			"{time=1, a=0},"
			"{time=2, a=0}]";

	if (this->is_enum_valued) {
		boost::regex e("=0");
		list_ok = boost::regex_replace(list_ok, e, "=\"white\"");
		list_ko = boost::regex_replace(list_ko, e, "=\"white\"");
	}

	this->field_.name("a");
	this->field_.set_input_list( this->input_list(list_ok) );

	this->field_.name("b");
	this->field_.set_input_list( this->input_list(list_ok) );

	this->field_.name("a");
	EXPECT_THROW_WHAT( {	this->field_.set_input_list( this->input_list(list_ko) );}, FieldCommonBase::ExcNonascendingTime, "for field 'a'" );
}


// check that it correctly introduce requered marks
TYPED_TEST(FieldFix, mark_input_times) {
	string list_ok = "["
			"{time=2, a=0}, {time=1, b=0}, {time=10, c=0},"
			"{time=3, a=0}, {time=2, b=0}, {time=1, c=0},"
			"{time=4, a=0}, {time=5, a=0, b=0}]";

	if (this->is_enum_valued) {
		boost::regex e("=0");
		list_ok = boost::regex_replace(list_ok, e, "=\"white\"");
	}

	this->field_.name("b");
	this->field_.set_input_list(this->input_list(list_ok));

	TimeGovernor tg;
	TimeMark::Type mark_type = tg.marks().new_mark_type();
	this->field_.mark_input_times(mark_type);
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
	EXPECT_ASSERT_DEATH( {this->set_field(this->my_domain, this->my_field_base);}, "Null; mesh pointer");

	EXPECT_EQ(nullptr, this->shared_->mesh_);
	this->set_mesh(*(this->my_mesh));
	EXPECT_EQ(this->my_mesh, this->shared_->mesh_);
	int size = this->my_mesh->region_db().size();
	EXPECT_EQ( size, this->region_fields_.size());
	EXPECT_EQ( size, this->data_->region_history_.size());


}


// Check copy constructor and assignment oprerator
TYPED_TEST(FieldFix, constructors) {
	typename TestFixture::FieldType f2(this->field_);
	typename TestFixture::FieldType f3;
	f3 = this->field_;
}

/*



TEST_F(FiledFix, get_const_accessor) {

}

*/

TYPED_TEST(FieldFix, set_field) {
	this->set_mesh(*(this->my_mesh));
	this->set_field(this->my_domain, this->my_field_base);

	Region reg = this->my_domain[0];
	auto const &history = this->data_->region_history_[reg.idx()];

	EXPECT_EQ(1, history.size());
	EXPECT_EQ(0.0, history[0].first);
	EXPECT_TRUE( bool(history[0].second) );

	this->set_field(this->my_domain, this->my_field_base, 3.0);
	EXPECT_EQ(2, history.size());
	EXPECT_EQ(3.0, history[0].first);
	EXPECT_EQ(0.0, history[1].first);

	EXPECT_ASSERT_DEATH( {this->set_field(this->my_domain, this->my_field_base, 1.0);}, "" );

	this->set_field(this->my_domain, this->my_field_base, 6.0);
	EXPECT_EQ(3, history.size());
	EXPECT_EQ(6.0, history[0].first);
	EXPECT_EQ(3.0, history[1].first);
	EXPECT_EQ(0.0, history[2].first);

	this->set_field(this->my_domain, this->my_field_base, 7.0);
	EXPECT_EQ(3, history.size());
	EXPECT_EQ(7.0, history[0].first);
	EXPECT_EQ(6.0, history[1].first);
	EXPECT_EQ(3.0, history[2].first);

}


TYPED_TEST(FieldFix, update_history) {
	string list_ok = "["
			"{time=0, r_set=\"ALL\", a =0, b =0},"
			"{time=1, r_set=\"BULK\", a =1, b =0},"
			"{time=2, r_set=\"BOUNDARY\", a =1, b =0},"
			"{time=3, r_set=\"ALL\", b =0},"
			"{time=4, r_set=\"ALL\", a =0},"
			"{time=5, r_set=\"ALL\", a =1}"
			"]";
	if (this->is_enum_valued) {
		DBGMSG("Enum Valued\n");
		list_ok = boost::regex_replace(list_ok, boost::regex(" =1"), "=\"white\"");
		list_ok = boost::regex_replace(list_ok, boost::regex(" =0"), "=\"black\"");
	}

	this->name("a");
	this->set_mesh(*(this->my_mesh));
	this->set_input_list( this->input_list(list_ok) );

	// time = 0.0
	TimeGovernor tg(0.0, 1.0);
	this->update_history(tg);

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

	double dt=tg.estimate_dt();
	tg.next_time();
	this->update_history(tg);

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
	this->update_history(tg);

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
	this->update_history(tg);

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
	this->update_history(tg);

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
	this->update_history(tg);

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


/*
TEST_F(FiledFix, set_time) {

}




TEST_F(FiledFix, check_initialized_region_field) {

}

*/
/*
namespace it = Input::Type;
TEST(Field, init_from_input) {
//    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
	Profiler::initialize();

	Mesh mesh;
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
	ifstream in(string( FilePath("mesh/simplest_cube.msh", FilePath::input_file) ).c_str());
	mesh.read_gmsh_from_stream(in);

    it::Selection sorption_type_sel =
            it::Selection("SorptionType")
            .add_value(1,"linear")
            .add_value(0,"none");


    Field<3, FieldValue<3>::Enum > sorption_type;
    Field<3, FieldValue<3>::Vector > init_conc;
    Field<3, FieldValue<3>::TensorFixed > conductivity;


    sorption_type.input_selection(sorption_type_sel);
    init_conc.n_comp(3);

    it::Record main_record =
            it::Record("main", "desc")
            .declare_key("sorption_type", sorption_type.make_input_tree(), it::Default::obligatory(), "desc")
            .declare_key("init_conc", init_conc.get_input_type(), it::Default::obligatory(), "desc")
            .declare_key("conductivity", conductivity.get_input_type(), it::Default::obligatory(), "desc");


    // read input string
    std::stringstream ss(field_input);
    Input::JSONToStorage reader;
    reader.read_stream( ss, main_record );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    sorption_type.set_mesh(mesh);
    init_conc.set_mesh(mesh);
    conductivity.set_mesh(mesh);

    auto r_set = mesh.region_db().get_region_set("BULK");

    sorption_type.set_field(r_set, in_rec.val<Input::AbstractRecord>("sorption_type"));
    init_conc.set_field(r_set, in_rec.val<Input::AbstractRecord>("init_conc"));
    conductivity.set_field(r_set, in_rec.val<Input::AbstractRecord>("conductivity"));


    sorption_type.set_time();
    init_conc.set_time();
    conductivity.set_time();

    {	

	    auto ele = mesh.element_accessor(5);

	    EXPECT_EQ( 1, sorption_type.value(ele.centre(), ele) );

	    EXPECT_TRUE( arma::min( arma::vec("10 20 30") == init_conc.value(ele.centre(), ele) ) );

	    arma::mat diff = arma::mat33("-0.5 0 0;0 0 0; 0 0 -0.5") - conductivity.value(ele.centre(), ele);
	    double norm=arma::norm(diff, 1);
	    EXPECT_DOUBLE_EQ( 0.0, norm );
    }

    {
	//  using const accessor
    	cout << "Second ele" << endl;
    	ElementAccessor<3> ele;

	    EXPECT_ASSERT_DEATH( {sorption_type.value(ele.centre(), ele);}  , "Invalid element accessor.");
	    Region reg = mesh.region_db().find_id(40);
	    EXPECT_TRUE( sorption_type.get_const_accessor(reg, ele));
	    EXPECT_TRUE( init_conc.get_const_accessor(reg, ele));
	    EXPECT_EQ( 1, sorption_type.value(ele.centre(), ele) );
	    EXPECT_TRUE( arma::min( arma::vec("10 20 30") == init_conc.value(ele.centre(), ele) ) );

   }



}








TEST(Field, init_from_default) {
//    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Profiler::initialize();
    
    Mesh mesh;
    ifstream in(string( FilePath("mesh/simplest_cube.msh", FilePath::input_file) ).c_str());
    mesh.read_gmsh_from_stream(in);

    Space<3>::Point p("1 2 3");

    {
        Field<3, FieldValue<3>::Scalar > scalar_field;

        // test default initialization of scalar field
        scalar_field.input_default( "45" );
        scalar_field.set_mesh(mesh);
        scalar_field.set_time();

        EXPECT_EQ( 45.0, scalar_field.value(p, mesh.element_accessor(0)) );
        EXPECT_EQ( 45.0, scalar_field.value(p, mesh.element_accessor(6)) );
        // this fails on dev.nti.tul.cz
        //EXPECT_DEATH( { scalar_field.value(p, mesh.element_accessor(0,true)); }, "Null field ptr " );
    }

    {
        Field<3, FieldValue<3>::Scalar > scalar_field("some", true);

        // test death of set_time without default value
        scalar_field.set_mesh(mesh);
        EXPECT_THROW_WHAT( {scalar_field.set_time();} , ExcXprintfMsg, "Missing value of the field");
    }
    //
    {
        Field<3, FieldValue<3>::Enum > enum_field("any", true);
        Input::Type::Selection sel("TestType");
        sel.add_value(0, "none")
           .add_value(1,"dirichlet")
           .close();

        enum_field.input_selection(sel);
        enum_field.input_default( "none" );
        enum_field.set_mesh(mesh);
        enum_field.set_time();

        EXPECT_EQ( 0 , enum_field.value(p, mesh.element_accessor(0, true)) );

    }
    Field<3, FieldValue<3>::Vector > vector_field;

}

/// Test optional fields dependent e.g. on BC type
TEST(Field, disable_where) {
//    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

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

    Mesh mesh;
    ifstream in(string( FilePath("mesh/simplest_cube.msh", FilePath::input_file) ).c_str());
    mesh.read_gmsh_from_stream(in);

    bc_type.set_mesh(mesh);
    bc_flux.set_mesh(mesh);
    bc_value.set_mesh(mesh);
    bc_sigma.set_mesh(mesh);

    /*
    1       37      "1D diagonal"
    2       38      "2D XY diagonal"
    2       101     ".top side"
    2       102     ".bottom side"
    3       39      "3D back"
    3       40      "3D front"
     */
/*

    typedef FieldConstant<3, FieldValue<3>::Scalar > SConst;
    typedef FieldConstant<3, FieldValue<3>::Enum > EConst;
    auto neumann_type = std::make_shared<EConst>();
    neumann_type->set_value(neumann);
    auto robin_type = std::make_shared<EConst>();
    robin_type->set_value(robin);
    auto one = std::make_shared<SConst>();
    one->set_value(1.0);

    bc_type.set_field(RegionSet(1, mesh.region_db().find_id(101)), neumann_type );
    bc_flux.set_field(RegionSet(1, mesh.region_db().find_id(101)), one );

    bc_type.set_field(RegionSet(1, mesh.region_db().find_id(102)), robin_type );
    bc_value.set_field(RegionSet(1, mesh.region_db().find_id(102)), one );
    bc_sigma.set_field(RegionSet(1, mesh.region_db().find_id(102)), one );

    bc_type.set_field(RegionSet(1, mesh.region_db().find_id(-3)), neumann_type );
    bc_flux.set_field(RegionSet(1, mesh.region_db().find_id(-3)), one );

    bc_type.set_time();
    bc_flux.set_time();
    bc_value.set_time();
    bc_sigma.set_time();
}


*/
