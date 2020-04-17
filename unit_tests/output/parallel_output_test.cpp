/*
 * parallel_output_test.cpp
 *
 *  Created on: Dec 7, 2018
 *      Author: David Flanderka
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include "config.h"

#include "mesh/mesh.h"
#include "io/output_time.hh"
#include "io/output_vtk.hh"
#include "io/output_mesh.hh"
#include "io/msh_gmshreader.h"
#include "system/logger_options.hh"
#include "system/sys_profiler.hh"
#include "fields/field_set.hh"
#include "fields/field.hh"
#include "la/distribution.hh"
#include "coupling/equation.hh"

#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"

namespace IT=Input::Type;


string eq_data_input = R"YAML(
data:
  - region: BULK
    time: 0.0
    init_scalar: !FieldFormula
      value: x
    init_vector: !FieldFormula
      value: 
        - x*y
        - x+y
        - z
)YAML";

/**
 * Special case of output test.
 *
 * It must run for n_proc=2 and tests collective of parallel computing result to serial output for continuous and discontinuous output.
 */
class TestParallelOutput : public testing::Test, public EquationBase
{
protected:
    class EqData : public FieldSet {
    public:
        EqData()
        {
            *this += init_scalar.name("init_scalar").description("Initial condition for scalar.").input_default("0.0");
            *this += init_vector.name("init_vector").description("Initial condition for vector.").input_default("0.0");
            init_scalar.units( UnitSI::dimensionless() );
            init_vector.units( UnitSI::dimensionless() );
        }

        Field<3, FieldValue<3>::Scalar > init_scalar;
        Field<3, FieldValue<3>::VectorFixed > init_vector;
    };

    class OutputVTKTest : public OutputVTK {
    public:
        OutputVTKTest() : OutputVTK() {};

        void gather_data(Mesh* mesh) {
            auto &elm_data_map = this->output_data_vec_[ELEM_DATA];
            for(unsigned int i=0; i<elm_data_map.size(); ++i) {
                auto serial_data = elm_data_map[i]->gather(mesh->get_el_ds(), mesh->get_el_4_loc());
                if (this->rank_==0) elm_data_map[i] = serial_data;
            }
        }

        std::shared_ptr<ElementDataCache<double>> nodes_cache() { return nodes_; };
        std::shared_ptr<ElementDataCache<unsigned int>> connectivity_cache() { return connectivity_; };
        std::shared_ptr<ElementDataCache<unsigned int>> offsets_cache() { return offsets_; };
        std::vector<OutputDataPtr> output_data_vec(OutputTime::DiscreteSpace space_type) { return this->output_data_vec_[space_type]; };
    };

    virtual void SetUp()
    {
        Profiler::instance();
        FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_2_elem.msh", FilePath::input_file);
        my_mesh = mesh_full_constructor("{mesh_file=\"" + (string)mesh_file + "\"}");

        component_names = { "comp_0", "comp_1", "comp_2" };
        stream = std::make_shared<OutputVTKTest>();
        rank = my_mesh->get_el_ds()->myp();
	}
    virtual void TearDown() {
        delete my_mesh;
    }

    IT::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( TestParallelOutput::EqData().make_field_descriptor_type("SomeEquation") )
                        .declare_key("init_scalar", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
                        .declare_key("init_vector", FieldAlgorithmBase< 3, FieldValue<3>::VectorFixed >::get_input_type_instance(), "" )
                        .close()
                        ), IT::Default::obligatory(), ""  )
                .close();
    }

    void check_distributions() {
        std::vector<int> expected_n_elems = {1, 1}; // expected numbers of elements for individual processes
        std::vector<int> expected_n_nodes = {3, 1}; // expected numbers of nodes for individual processes
        EXPECT_EQ(this->my_mesh->get_el_ds()->lsize(), expected_n_elems[rank]);
        EXPECT_EQ(this->my_mesh->get_node_ds()->lsize(), expected_n_nodes[rank]);
    }

    void read_input(const string &input) {
        // read input string
        Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_YAML );
        Input::Record in_rec=reader.get_root_interface<Input::Record>();

        TimeGovernor tg(0.0, 1.0);

        data.set_components(component_names);        // set number of substances posibly read from elsewhere

        DebugOut() << "init\n";

        static std::vector<Input::Array> inputs;
        unsigned int input_last = inputs.size(); // position of new item
        inputs.push_back( in_rec.val<Input::Array>("data") );

        data.set_mesh(*my_mesh);
        data.set_input_list( inputs[input_last], tg );
        data.set_time(tg.step(), LimitSide::right);
    }

    void make_output_mesh() {
        output_mesh->create_sub_mesh();
        output_mesh->make_serial_master_mesh();
        stream->set_output_data_caches(output_mesh);
    }

	Mesh * my_mesh;
	std::shared_ptr<OutputMeshBase> output_mesh;
    EqData data;
    std::vector<string> component_names;
    std::shared_ptr<OutputVTKTest> stream;
    int rank;
};



TEST_F(TestParallelOutput, continuous_mesh)
{
    check_distributions();
    read_input(eq_data_input);
    output_mesh = std::make_shared<OutputMesh>(*my_mesh);
    make_output_mesh();

    /* Tests of output data */
    if (rank == 0) {
        std::vector<double> expected_nodes = { -3, 3, 0, 3, 3, 0, 3, -3, 0, -3, -3, 0 };
        std::vector<unsigned int> expected_connectivities = { 0, 3, 2, 0, 2, 1 };
        std::vector<unsigned int> expected_offsets = { 3, 6 };
        auto &node_vec = *( stream->nodes_cache()->get_component_data(0).get() );
        for (unsigned int i=0; i<node_vec.size(); ++i) EXPECT_DOUBLE_EQ( node_vec[i], expected_nodes[i]);
        auto &conn_vec = *( stream->connectivity_cache()->get_component_data(0).get() );
        for (unsigned int i=0; i<conn_vec.size(); ++i) EXPECT_EQ( conn_vec[i], expected_connectivities[i]);
        auto &offs_vec = *( stream->offsets_cache()->get_component_data(0).get() );
        for (unsigned int i=0; i<offs_vec.size(); ++i) EXPECT_EQ( offs_vec[i], expected_offsets[i]);

        std::map<string, std::vector<double>> expected_field_data;
        expected_field_data["init_scalar"] = { -1, 1 };
        expected_field_data["init_vector"] = { 1, -2, 0, 1, 2, 0 };
        auto elem_cache_map = stream->output_data_vec(OutputTime::DiscreteSpace::ELEM_DATA);
        for (unsigned int i=0; i<elem_cache_map.size(); ++i) {
            std::shared_ptr< ElementDataCache<double> > current_cache = dynamic_pointer_cast<ElementDataCache<double> >(elem_cache_map[i]);
            string field_name = current_cache->field_input_name();
            auto &current_vec = *( current_cache->get_component_data(0).get() );
            EXPECT_EQ( expected_field_data[field_name].size(), current_vec.size() );
            for (unsigned int i=0; i<current_vec.size(); ++i) EXPECT_DOUBLE_EQ( expected_field_data[field_name][i], current_vec[i] );
        }
    }
}


TEST_F(TestParallelOutput, discontinuous_mesh)
{
    check_distributions();
    read_input(eq_data_input);
    output_mesh = std::make_shared<OutputMeshDiscontinuous>(*my_mesh);
    make_output_mesh();

    /* Simulate field output */
    data.init_scalar.field_output(stream);
    data.init_vector.field_output(stream);
    stream->gather_data(my_mesh);

    /* Tests of output data */
    if (rank == 0) {
        std::vector<double> expected_nodes = { -3, 3, 0, -3, -3, 0, 3, -3, 0, -3, 3, 0, 3, -3, 0, 3, 3, 0 };
        std::vector<unsigned int> expected_connectivities = { 0, 1, 2, 3, 4, 5 };
        std::vector<unsigned int> expected_offsets = { 3, 6 };
        auto &node_vec = *( stream->nodes_cache()->get_component_data(0).get() );
        for (unsigned int i=0; i<node_vec.size(); ++i) EXPECT_DOUBLE_EQ( node_vec[i], expected_nodes[i]);
        auto &conn_vec = *( stream->connectivity_cache()->get_component_data(0).get() );
        for (unsigned int i=0; i<conn_vec.size(); ++i) EXPECT_EQ( conn_vec[i], expected_connectivities[i]);
        auto &offs_vec = *( stream->offsets_cache()->get_component_data(0).get() );
        for (unsigned int i=0; i<offs_vec.size(); ++i) EXPECT_EQ( offs_vec[i], expected_offsets[i]);

        std::map<string, std::vector<double>> expected_field_data;
        expected_field_data["init_scalar"] = { -1, 1 };
        expected_field_data["init_vector"] = { 1, -2, 0, 1, 2, 0 };
        auto elem_cache_map = stream->output_data_vec(OutputTime::DiscreteSpace::ELEM_DATA);
        for (unsigned int i=0; i<elem_cache_map.size(); ++i) {
            std::shared_ptr< ElementDataCache<double> > current_cache = dynamic_pointer_cast<ElementDataCache<double> >(elem_cache_map[i]);
            string field_name = current_cache->field_input_name();
            auto &current_vec = *( current_cache->get_component_data(0).get() );
            EXPECT_EQ( expected_field_data[field_name].size(), current_vec.size() );
            for (unsigned int i=0; i<current_vec.size(); ++i) EXPECT_DOUBLE_EQ( expected_field_data[field_name][i], current_vec[i] );
        }
    }
}
