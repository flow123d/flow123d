/*
 * observe_test.cpp
 *
 *  Created on: Jun 29, 2016
 *      Author: jb
 */



#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include "flow_gtest_mpi.hh"
#include "io/observe.hh"
#include "mesh/mesh.h"
#include "input/reader_to_storage.hh"
#include "input/accessors.hh"
#include "system/sys_profiler.hh"
#include "fields/field_set.hh"
#include "armadillo"


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
    { point: [0, -0.5, -0.5] },    // element 7
    { name: "s_el7",point: [0, -0.5, -0.5], snap_dim: 3 },
    { name: "s_front_face_el7", point: [0, -1, -0.5], snap_dim: 2 },
    { name: "s_front_bot_edge_el7", point: [0.2, -0.9, -0.9], snap_dim: 1 },
    { name: "s_2d_el2", point: [0, -0.5, -0.5], snap_region: "2D XY diagonal" },
    { name: "s_2d_el2", point: [0, -0.5, -0.5], snap_region: "2D XY diagonal", snap_dim: 2}
  ],
  observe_fields: [
  ]
}
)JSON";


class TestObserve;

class TestObservePoint : public ObservePoint {
    TestObservePoint(const ObservePoint &point)
    : ObservePoint(point)
    {}

    friend TestObserve;
};

class TestObserve : public Observe {
public:
    TestObserve(Mesh &mesh, Input::Record in_rec)
    : Observe("test_eq", mesh, in_rec)
    {}

    void check_points_input() {
        vector<TestObservePoint> my_points;
        for(auto &point: this->points_) my_points.push_back(TestObservePoint(point));

        EXPECT_EQ("obs_0", my_points[0].name_);


        auto vec_match = arma::vec3({0.0, -0.5, -0.5}) == my_points[0].input_point_;
        EXPECT_EQ(1, arma::min( vec_match) );
        EXPECT_EQ(4, my_points[0].snap_dim_);
        EXPECT_EQ("ALL", my_points[0].snap_region_);

        EXPECT_EQ("s_el7", my_points[1].name_);
        EXPECT_EQ(3, my_points[1].snap_dim_);

        EXPECT_EQ(2, my_points[2].snap_dim_);
        EXPECT_EQ("s_front_face_el7", my_points[2].name_);
    }

};



class EqData : public FieldSet {
public:
    typedef Field<3, FieldValue<3>::Scalar > ScalarField;
    typedef Field<3, FieldValue<3>::Enum > EnumField;
    typedef Field<3, FieldValue<3>::VectorFixed > VectorField;
    typedef Field<3, FieldValue<2>::TensorFixed > TensorField;

    EqData() {
        ADD_FIELD(scalar_field, "").units(UnitSI::one());
        ADD_FIELD(enum_field, "").units(UnitSI::one());
        ADD_FIELD(vector_field, "").units(UnitSI::one());
        ADD_FIELD(tensor_field, "").units(UnitSI::one());
    }

    ScalarField scalar_field;
    EnumField enum_field;
    VectorField vector_field;
    TensorField tensor_field;
};



TEST(Observe, all) {
    Profiler::initialize();
    EqData field_set;

    auto field_selection = field_set.make_output_field_selection("eq_data", "").close();
    auto output_type = Input::Type::Record("Output", "")
        .declare_key("observe_fields", Input::Type::Array(field_selection), Input::Type::Default::obligatory(), "" )
        .declare_key("observe_points", Input::Type::Array(ObservePoint::get_input_type()), Input::Type::Default::obligatory(), "")
        .close();
    auto in_rec = Input::ReaderToStorage(test_input, output_type, Input::FileFormat::format_JSON)
        .get_root_interface<Input::Record>();

    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh", FilePath::input_file);
    Mesh *mesh = new Mesh();
    ifstream in(string(mesh_file).c_str());
    mesh->read_gmsh_from_stream(in);

    TestObserve obs(*mesh, in_rec);
    obs.check_points_input();


}
