/*
 * observe_test.cpp
 *
 *  Created on: Jun 29, 2016
 *      Author: jb
 */



#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS



#include "flow_gtest_mpi.hh"
#include <mesh_constructor.hh>
#include "io/observe.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "input/reader_to_storage.hh"
#include "input/accessors.hh"
#include "system/sys_profiler.hh"
#include "fields/field_set.hh"
#include "fields/field.hh"
#include "fields/equation_output.hh"
#include "fields/assembly_observe.hh"
#include "armadillo"
#include "system/armadillo_tools.hh"
#include "tools/time_governor.hh"
#include "fem/dofhandler.hh"
#include "fem/discrete_space.hh"
#include "fem/fe_p.hh"
#include "../arma_expect.hh"
#include <fstream>




FLOW123D_FORCE_LINK_IN_PARENT(field_constant)
FLOW123D_FORCE_LINK_IN_PARENT(field_formula)

/**
 * Supposed to be used on 'simplest_cube':
 * box [-1,-1,-1] [1,1,1]
 *  right face X=1, front face is Y=-1, top face is Z=1

$PhysicalNames
6
1       37      "1D diagonal"
2       38      "2D XY diagonal"
2       101     ".top side"
2       102     ".bottom side"
3       39      "3D back"
3       40      "3D front"
$EndPhysicalNames
 */

const string test_input = R"JSON(
{
  observe_points: [
    // observe element 7, various snapping
    { point: [0, -0.5, -0.5] },    
    { name: "s_el7",point: [0, -0.6, -0.6], snap_dim: 3 },
    { name: "s_front_face_el7", point: [0, -1, -0.5], snap_dim: 2 },
    { name: "s_front_bot_edge_el7", point: [0.2, -0.9, -0.9], snap_dim: 1 },
    { name: "s_front_node_el7", point: [-0.8, -0.9, -0.9], snap_dim: 0 },
    // find 2D observe element 2, various snapping
    { name: "s_2d_el2", point: [0, -0.5, -0.5], snap_region: "2D XY diagonal" },
    { name: "s_2d_el2", point: [0, -0.5, -0.5], snap_region: "2D XY diagonal", snap_dim: 2},
    // find 1D observe element , snap to node [-1, -1 ,1]
    { name: "s_1d_el2", point: [-0.5, -0.5, 0], snap_region: "1D diagonal", snap_dim: 0}
  ],
  input_fields: [
       {
         region: "ALL",
         scalar_field: {TYPE: "FieldFormula", value: "x+y+2*z"}, 
         enum_field: "one",
         vector_field: [0,2,3],
         tensor_field: [1, 2, 3, 4, 5, 6]
       }
  ]
}
)JSON";


class TestObserve;

class TestObservePoint : public ObservePoint {
public:
    TestObservePoint(const ObservePoint &point)
    : ObservePoint(point)
    {}

    // Constructor simulate input from YAML file.
    TestObservePoint(string point_str, unsigned int snap, string region_name)
    : ObservePoint()
    {
        observe_data_.distance_ = numeric_limits<double>::infinity();
        input_point_= arma::vec3(point_str);
        snap_dim_ = snap;
        snap_region_name_ = region_name;
        max_search_radius_ = 0.1;
    }

    // Check find_observe_point.
    void check(Mesh &mesh, string local_str, string global_point_str, unsigned int i_elm) {
        find_observe_point(mesh);
        EXPECT_EQ(i_elm, observe_data_.element_idx_);
        if (local_str != "") {
            EXPECT_ARMA_EQ( arma::vec(local_str), observe_data_.local_coords_);
        }
        EXPECT_ARMA_EQ( arma::vec3(global_point_str), observe_data_.global_coords_);
    }

    friend TestObserve;
};

class TestObserve : public Observe {
public:
    TestObserve(Mesh &mesh, Input::Array in_array)
    : Observe("test_eq", mesh, in_array, 6, std::make_shared<TimeUnitConversion>())
    {
        for(auto &point: this->points_) my_points.push_back(TestObservePoint(point));
    }

    void check_points_input() {

        EXPECT_EQ("obs_0", my_points[0].name_);

        EXPECT_ARMA_EQ(arma::vec3({0.0, -0.5, -0.5}), my_points[0].input_point_ );
        EXPECT_EQ(4, my_points[0].snap_dim_);
        EXPECT_EQ("ALL", my_points[0].snap_region_name_);

        EXPECT_EQ("s_el7", my_points[1].name_);
        EXPECT_EQ(3, my_points[1].snap_dim_);

        EXPECT_EQ(2, my_points[2].snap_dim_);
        EXPECT_EQ("s_front_face_el7", my_points[2].name_);
    }

    void check_observe_points() {
        // no snap, [0, -0.5, -0.5]
        EXPECT_EQ(6,  my_points[0].observe_data_.element_idx_);
        EXPECT_ARMA_EQ( arma::vec3("0 -0.5 -0.5"),  my_points[0].observe_data_.global_coords_);
        EXPECT_ARMA_EQ( arma::vec("0.25 0.25 0.25"),  my_points[0].observe_data_.local_coords_);
        EXPECT_DOUBLE_EQ(0, my_points[0].observe_data_.distance_);

        // snap 3, [0, -0.6, -0.6]
        EXPECT_EQ(6,  my_points[1].observe_data_.element_idx_);
        EXPECT_ARMA_EQ( arma::vec3("0 -0.5 -0.5"),  my_points[1].observe_data_.global_coords_);
        EXPECT_ARMA_EQ( arma::vec("0.25 0.25 0.25"),  my_points[1].observe_data_.local_coords_);


        // snap 2, [0, -1, -0.5]
        EXPECT_EQ(6,  my_points[2].observe_data_.element_idx_);
        EXPECT_ARMA_EQ( arma::vec3({-1.0/3, -1, -1.0/3}),  my_points[2].observe_data_.global_coords_);
        EXPECT_ARMA_EQ( arma::vec({1.0/3, 1.0/3, 0}),  my_points[2].observe_data_.local_coords_);

        // snap 1, [0.2, -0.9, -0.9]
        EXPECT_EQ(6,  my_points[3].observe_data_.element_idx_);
        EXPECT_ARMA_EQ( arma::vec3({0, -1, -1}),  my_points[3].observe_data_.global_coords_);
        EXPECT_ARMA_EQ( arma::vec({0.5, 0.5, 0}),  my_points[3].observe_data_.local_coords_);

        // snap 0, [-0.8, -0.9, -0.9]
        EXPECT_EQ(6,  my_points[4].observe_data_.element_idx_);
        EXPECT_ARMA_EQ( arma::vec3({-1, -1, -1}),  my_points[4].observe_data_.global_coords_);
        EXPECT_ARMA_EQ( arma::vec({0, 1, 0}),  my_points[4].observe_data_.local_coords_);

        // find 2D observe element 2, various snapping
        //{ name: "s_2d_el2", point: [0, -0.5, -0.5], snap_region: "2D XY diagonal" },
        EXPECT_EQ(1,  my_points[5].observe_data_.element_idx_);
        EXPECT_ARMA_EQ( arma::vec3({-1.0/4, -1.0/4, -1.0/2}),  my_points[5].observe_data_.global_coords_);
        EXPECT_ARMA_EQ( arma::vec({3.0/8, 3.0/8}),  my_points[5].observe_data_.local_coords_);


        //{ name: "s_2d_el2", point: [0, -0.5, -0.5], snap_region: "2D XY diagonal", snap_dim: 2},
        EXPECT_EQ(1,  my_points[6].observe_data_.element_idx_);
        EXPECT_ARMA_EQ( arma::vec3({-1.0/3, -1.0/3, -1.0/3}),  my_points[6].observe_data_.global_coords_);
        EXPECT_ARMA_EQ( arma::vec({1.0/3, 1.0/3}),  my_points[6].observe_data_.local_coords_);

        // find 1D observe element , snap to node [-1, -1 ,1]
        //{ name: "s_1d_el1", point: [-0.5, -0.5, 0], snap_region: "1D diagonal", snap_dim: 0}
        EXPECT_EQ(0,  my_points[7].observe_data_.element_idx_);
        EXPECT_ARMA_EQ( arma::vec3({-1, -1, 1}),  my_points[7].observe_data_.global_coords_);
        EXPECT_ARMA_EQ( arma::vec({0}),  my_points[7].observe_data_.local_coords_);


    }
    std::vector<TestObservePoint> my_points;

};



class EqData : public EquationOutput {
public:
    typedef Field<3, FieldValue<3>::Scalar > ScalarField;
    typedef Field<3, FieldValue<3>::Enum > EnumField;
    typedef Field<3, FieldValue<3>::VectorFixed > VectorField;
    typedef Field<3, FieldValue<3>::TensorFixed > TensorField;

    EqData() {
        static Input::Type::Selection  selection =
                Input::Type::Selection("test_enum")
                     .add_value(0, "zero")
                     .add_value(1, "one")
                     .close();

    
        *this += scalar_field.name("scalar_field").units(UnitSI::one());
        *this += vector_field.name("vector_field").units(UnitSI::one());
        *this += tensor_field.name("tensor_field").units(UnitSI::one());
        *this += enum_field.name("enum_field").units(UnitSI::one()).input_selection(selection);

        this->add_coords_field();
        this->set_default_fieldset();
    }

    ScalarField scalar_field;
    EnumField enum_field;
    VectorField vector_field;
    TensorField tensor_field;
};


TEST(ObservePoint, find_observe_point) {
    Profiler::instance();
    armadillo_setup();

    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh", FilePath::input_file);
    Mesh *mesh = mesh_full_constructor("{ mesh_file=\"" + (string)mesh_file + "\", optimize_mesh=false }");

    //DebugOut() << "obs1";
    auto obs = TestObservePoint("0 -0.5 -0.5", 4, "ALL");
    obs.check(*mesh,"0.25 0.25 0.25", "0 -0.5 -0.5", 6);

    //DebugOut() << "obs2";
    auto obs2 = TestObservePoint("1 1 1.001", 4, "3D front");
    obs2.check(*mesh,"0 0 0", "1 1 1", 8);
    Profiler::uninitialize();
}


TEST(Observe, all) {
    Profiler::instance();
    armadillo_setup();
    std::shared_ptr<EqData> field_set = std::make_shared<EqData>();

    auto output_type = Input::Type::Record("Output", "")
        .declare_key("observe_points", Input::Type::Array(ObservePoint::get_input_type()), Input::Type::Default::obligatory(), "")
        .declare_key("input_fields", Input::Type::Array(
                EqData()
                .make_field_descriptor_type("SomeEquation")
                .close() ), Input::Type::Default::obligatory(), "")
        .close();
    auto in_rec = Input::ReaderToStorage(test_input, output_type, Input::FileFormat::format_JSON)
        .get_root_interface<Input::Record>();

    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh", FilePath::input_file);
    Mesh *mesh = mesh_full_constructor("{ mesh_file=\"" + (string)mesh_file + "\", optimize_mesh=false, global_snap_radius=1.0 }");
    field_set->set_mesh(*mesh);
    MixedPtr<FE_P_disc> fe_p_disc(0);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe_p_disc);
    dh->distribute_dofs(ds);

    {
    std::shared_ptr<TestObserve> obs = std::make_shared<TestObserve>(*mesh, in_rec.val<Input::Array>("observe_points"));
    obs->check_points_input();
    obs->check_observe_points();

    std::unordered_set<string> observe_fields_list;
    observe_fields_list.insert("scalar_field");
    observe_fields_list.insert("vector_field");
    observe_fields_list.insert("tensor_field");
    //observe_fields_list.insert("enum_field");
    GenericAssemblyObserve< AssemblyObserveOutput > * observe_output_assembly;
    observe_output_assembly = new GenericAssemblyObserve< AssemblyObserveOutput >( field_set.get(), observe_fields_list, std::dynamic_pointer_cast<Observe>(obs) );

    // read fiels
    TimeGovernor tg(0.0, 1.0);
    field_set->set_input_list( in_rec.val<Input::Array>("input_fields"), tg );
    field_set->set_time(tg.step(), LimitSide::right);

    for(auto field_name : observe_fields_list) {
        auto &field = (*field_set)[field_name];
        obs->prepare_compute_data(field.name(), field.time(), field.n_shape());
    }

    observe_output_assembly->assemble(dh);
    obs->output_time_frame( true );

    tg.next_time();
    field_set->set_time( tg.step(), LimitSide::right);
    observe_output_assembly->assemble(dh);
    }
    // closed observe file 'test_eq_observe.yaml'
    // check results
    std::ifstream  obs_file("test_eq_observe.yaml");
    std::stringstream str_obs_file;
    str_obs_file << obs_file.rdbuf();
    obs_file.close();

    std::ifstream  obs_file_ref("test_eq_observe_ref.yaml");
    std::stringstream str_obs_file_ref;
    str_obs_file_ref << obs_file_ref.rdbuf();
    obs_file_ref.close();

//    if (mesh->get_el_ds()->myp()==0)
//        EXPECT_EQ(str_obs_file_ref.str(), str_obs_file.str());
    Profiler::uninitialize();
}

