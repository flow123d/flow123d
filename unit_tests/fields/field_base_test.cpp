/*
 * field_base_test.cpp
 *
 *  Created on: Feb 3, 2013
 *      Author: jb
 */

#define TEST_USE_MPI

#include <memory>
#include <flow_gtest_mpi.hh>

#include "fields/fied.hh"
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


    sorption_type.set_selection(&sorption_type_sel);
    init_conc.set_n_comp(3);

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
        scalar_field.set_default( "45" );
        scalar_field.set_mesh(mesh);
        scalar_field.set_time();

        EXPECT_EQ( 45.0, scalar_field.value(p, mesh.element_accessor(0)) );
        EXPECT_EQ( 45.0, scalar_field.value(p, mesh.element_accessor(6)) );
        // this fails on dev.nti.tul.cz
        //EXPECT_DEATH( { scalar_field.value(p, mesh.element_accessor(0,true)); }, "Null field ptr " );
    }

    {
        BCField<3, FieldValue<3>::Scalar > scalar_field;

        // test death of set_time without default value
        scalar_field.set_mesh(mesh);
        EXPECT_THROW_WHAT( {scalar_field.set_time();} , ExcXprintfMsg, "Missing value of the field");
    }
    //
    {
        BCField<3, FieldValue<3>::Enum > enum_field;
        Input::Type::Selection sel("TestType");
        sel.add_value(0, "none")
           .add_value(1,"dirichlet")
           .close();

        enum_field.set_selection(&sel);
        enum_field.set_default( "none" );
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
    BCField<3, FieldValue<3>::Enum > bc_type;
    bc_type.set_name("bc_type");

    std::vector<FieldEnum> list;
    BCField<3, FieldValue<3>::Scalar > bc_value;
    bc_value.set_name("bc_value");
    list.clear(); list.push_back(neumann);
    bc_value.disable_where( &bc_type, list );

    BCField<3, FieldValue<3>::Scalar > bc_flux;
    bc_flux.set_name("bc_flux");
    list.clear(); list.push_back(dirichlet); list.push_back(robin);
    bc_flux.disable_where( &bc_type, list );

    BCField<3, FieldValue<3>::Scalar > bc_sigma;
    bc_sigma.set_name("bc_sigma");
    list.clear(); list.push_back(dirichlet); list.push_back(neumann);
    bc_sigma.disable_where( &bc_type, list );

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

