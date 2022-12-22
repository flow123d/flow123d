/*
 * multi_field_test.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include <arma_expect.hh>

#include <fields/multi_field.hh>
#include <fields/field_fe.hh>
#include <fields/field_set.hh>
#include <fields/eval_points.hh>
#include <fields/eval_subset.hh>
#include <fields/field_value_cache.hh>
#include <fields/field_set.hh>
#include <fields/field_flag.hh>
#include <coupling/generic_assembly.hh>
#include <coupling/assembly_base.hh>
#include <quadrature/quadrature.hh>
#include <quadrature/quadrature_lib.hh>
#include <io/msh_gmshreader.h>
#include <fem/dofhandler.hh>
#include <fem/dh_cell_accessor.hh>
#include <mesh/mesh.h>
#include <mesh/accessors.hh>
#include <mesh/range_wrapper.hh>
#include <tools/unit_si.hh>
#include <input/type_base.hh>
#include <input/reader_to_storage.hh>
#include <input/type_output.hh>
#include <system/sys_profiler.hh>


#include <iostream>
using namespace std;

FLOW123D_FORCE_LINK_IN_PARENT(field_constant)
FLOW123D_FORCE_LINK_IN_PARENT(field_formula)
FLOW123D_FORCE_LINK_IN_PARENT(field_fe)

string field_constant_input = R"YAML(
common: !FieldConstant 
  value:
   - 1
   - 2
   - 3
transposed:
 - !FieldConstant
   value: 1
 - !FieldConstant
   value: 2
 - !FieldConstant
   value: 3
)YAML";

TEST(MultiField, transposition) {
	MultiField<3, FieldValue<3>::Scalar> empty_mf;
	Input::Type::Record in_rec = Input::Type::Record("MultiFieldTest","")
	    .declare_key("common", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
	    .declare_key("transposed", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
	    .close();

	Input::ReaderToStorage json_reader(field_constant_input, in_rec, Input::FileFormat::format_YAML);
	Input::Record input = json_reader.get_root_interface<Input::Record>();

	Input::Array common = input.val<Input::Array>("common");
	Input::Array transposed = input.val<Input::Array>("transposed");

	EXPECT_EQ(common.size(), transposed.size());

	auto it_c = common.begin<Input::AbstractRecord>();
	for (auto it_t = transposed.begin<Input::AbstractRecord>(); it_t != transposed.end(); ++it_t, ++it_c) {
		Input::Record rec_t = (*it_t);
		Input::Record rec_c = (*it_c);
		EXPECT_DOUBLE_EQ( rec_t.val<double>("value"), rec_c.val<double>("value") );
	}
}



class MultiFieldTest : public testing::Test {
public:

    class EqData : public FieldSet, public ElementCacheMap {
    public:
		EqData() {
            *this+=scalar_field
            		.name("scalar_field")
        			.description("Scalar field.")
        			.units(UnitSI().m())
					.input_default("0.0")
        			.flags_add( in_main_matrix );

            *this+=vector_field
            		.name("vector_field")
        			.description("Vector field.")
        			.units(UnitSI().m())
					.input_default("0.0")
        			.flags_add( in_main_matrix );

            // Asumme following types:
            eval_points_ = std::make_shared<EvalPoints>();
            Quadrature *q_bulk_1d = new QGauss(1, 0);
            Quadrature *q_bulk_2d = new QGauss(2, 0);
            Quadrature *q_bulk_3d = new QGauss(3, 0);
            bulk_int[0] = eval_points_->add_bulk<1>(*q_bulk_1d );
            bulk_int[1] = eval_points_->add_bulk<2>(*q_bulk_2d );
            bulk_int[2] = eval_points_->add_bulk<3>(*q_bulk_3d );
            this->init(eval_points_);
        }

        void register_eval_points() {
            unsigned int reg_idx = computed_dh_cell_.elm().region_idx().idx();
            for (auto p : bulk_int[computed_dh_cell_.dim()-1]->points(this->position_in_cache(computed_dh_cell_.elm_idx()), this) ) {
                this->eval_point_data_.emplace_back(reg_idx, computed_dh_cell_.elm_idx(), p.eval_point_idx(), computed_dh_cell_.local_idx());
            }
            this->eval_point_data_.make_permanent();
        }

        void update_cache() {
            this->register_eval_points();
            this->create_patch();
            this->cache_update(*this);
            this->finish_elements_update();
        }


        // fields
        MultiField<3, FieldValue<3>::Scalar > scalar_field;
        MultiField<3, FieldValue<3>::VectorFixed > vector_field;
        std::shared_ptr<EvalPoints> eval_points_;
        std::array<std::shared_ptr<BulkIntegral>, 3> bulk_int;  // dim 1,2,3
        std::shared_ptr<DOFHandlerMultiDim> dh_;
        DHCellAccessor computed_dh_cell_;
    };

    MultiFieldTest() : tg(0.25, 0.75) {
    	Profiler::instance();
    	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    	component_names = {"component_0", "component_1", "component_2"};

        eq_data_ = std::make_shared<EqData>();
        eq_data_->add_coords_field();
        eq_data_->set_default_fieldset();
        mesh_ = mesh_full_constructor("{ mesh_file=\"fields/simplest_cube_data.msh\", optimize_mesh=false }");
        eq_data_->dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

	~MultiFieldTest() {
    	delete mesh_;
    }

    void read_input(const string &input) {
        // read input string
        Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_YAML );
        Input::Record in_rec=reader.get_root_interface<Input::Record>();

        eq_data_->set_components(component_names);        // set number of substances possibly read from elsewhere

        static std::vector<Input::Array> inputs;
        unsigned int input_last = inputs.size(); // position of new item
        inputs.push_back( in_rec.val<Input::Array>("data") );

        eq_data_->set_mesh(*mesh_);
        eq_data_->set_input_list( inputs[input_last], tg );
    }

    void set_dh_cell(unsigned int elm_idx, unsigned int reg_idx) {
        ElementAccessor<3> elm = mesh_->element_accessor(elm_idx);
        EXPECT_EQ(reg_idx, elm.region().id());                         // check element accessor
        eq_data_->computed_dh_cell_ = this->eq_data_->dh_->cell_accessor_from_element(elm.idx());
    }


    static Input::Type::Record & get_input_type();
    static MultiField<3, FieldValue<3>::Scalar> empty_mf;
    std::shared_ptr<EqData> eq_data_;
    std::vector<std::string> component_names;
    Mesh * mesh_;
    TimeGovernor tg;
};

MultiField<3, FieldValue<3>::Scalar> MultiFieldTest::empty_mf = MultiField<3, FieldValue<3>::Scalar>();

Input::Type::Record & MultiFieldTest::get_input_type() {
    return IT::Record("SomeEquation","")
            .declare_key("data", IT::Array(
                IT::Record("MultiField_Data", FieldCommon::field_descriptor_record_description("MultiField_Data") )
                .copy_keys( MultiFieldTest::EqData().make_field_descriptor_type("MultiField") )
                .declare_key("scalar_field", empty_mf.get_multifield_input_type(), "" )
                .declare_key("vector_field", empty_mf.get_multifield_input_type(), "" )
                .close()
                ), IT::Default::obligatory(), ""  )
            .close();
}


TEST_F(MultiFieldTest, const_full_test) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldConstant
          value:
           - 1
           - 2
           - 3
    )YAML";

    this->read_input(eq_data_input);
    this->set_dh_cell(4, 39);

    for (uint i_time=0; i_time<2; i_time++) { // test in 2 time steps: 0.25, 1.0
        eq_data_->set_time(tg.step(), LimitSide::right);
        eq_data_->cache_reallocate( *(eq_data_.get()), *(eq_data_.get()) );
        eq_data_->update_cache();
        auto p = *( eq_data_->bulk_int[0]->points(eq_data_->position_in_cache(eq_data_->computed_dh_cell_.elm_idx()), eq_data_.get()).begin() );

        EXPECT_EQ(3, eq_data_->scalar_field.size());
        for (uint i_comp=0; i_comp<3; ++i_comp) {
            EXPECT_EQ(i_comp+1, eq_data_->scalar_field[i_comp](p));
        }

        tg.next_time();
    }
}


TEST_F(MultiFieldTest, const_base_test) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldConstant
          value: 1
    )YAML";

    this->read_input(eq_data_input);
    this->set_dh_cell(4, 39);

    for (uint i_time=0; i_time<2; i_time++) { // test in 2 time steps: 0.25, 1.0
        eq_data_->set_time(tg.step(), LimitSide::right);
        eq_data_->cache_reallocate( *(eq_data_.get()), *(eq_data_.get()) );
        eq_data_->update_cache();
        auto p = *( eq_data_->bulk_int[0]->points(eq_data_->position_in_cache(eq_data_->computed_dh_cell_.elm_idx()), eq_data_.get()).begin() );

        EXPECT_EQ(3, eq_data_->scalar_field.size());
        for (uint i_comp=0; i_comp<3; ++i_comp) {
            EXPECT_EQ(1, eq_data_->scalar_field[i_comp](p));
        }

        tg.next_time();
    }
}


TEST_F(MultiFieldTest, const_autoconv_test) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: 1
    )YAML";

    this->read_input(eq_data_input);
    this->set_dh_cell(4, 39);

    for (uint i_time=0; i_time<2; i_time++) { // test in 2 time steps: 0.25, 1.0
        eq_data_->set_time(tg.step(), LimitSide::right);
        eq_data_->cache_reallocate( *(eq_data_.get()), *(eq_data_.get()) );
        eq_data_->update_cache();
        auto p = *( eq_data_->bulk_int[0]->points(eq_data_->position_in_cache(eq_data_->computed_dh_cell_.elm_idx()), eq_data_.get()).begin() );

        EXPECT_EQ(3, eq_data_->scalar_field.size());
        for (uint i_comp=0; i_comp<3; ++i_comp) {
            EXPECT_EQ(1, eq_data_->scalar_field[i_comp](p));
        }

        tg.next_time();
    }
}


TEST_F(MultiFieldTest, formula_full_test) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFormula
          value:
           - t
           - x
           - y-t
    )YAML";

    this->read_input(eq_data_input);
    this->set_dh_cell(4, 39);

    for (uint i_time=0; i_time<2; i_time++) { // test in 2 time steps: 0.25, 1.0
        eq_data_->set_time(tg.step(), LimitSide::right);
        eq_data_->cache_reallocate( *(eq_data_.get()), *(eq_data_.get()) );
        eq_data_->update_cache();
        auto p = *( eq_data_->bulk_int[0]->points(eq_data_->position_in_cache(eq_data_->computed_dh_cell_.elm_idx()), eq_data_.get()).begin() );

        arma::vec3 elm_cntr = this->eq_data_->computed_dh_cell_.elm().centre(); // {-0.5, 0.5, 0.0}

        // test of FieldFormula - full input
        EXPECT_EQ(tg.t(), eq_data_->scalar_field[0](p));
        EXPECT_EQ(elm_cntr(0), eq_data_->scalar_field[1](p));
        EXPECT_EQ(elm_cntr(1)-tg.t(), eq_data_->scalar_field[2](p));

        tg.next_time();
    }
}


TEST_F(MultiFieldTest, formula_base_test) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFormula
          value: x
    )YAML";

    this->read_input(eq_data_input);
    this->set_dh_cell(4, 39);

    for (uint i_time=0; i_time<2; i_time++) { // test in 2 time steps: 0.25, 1.0
        eq_data_->set_time(tg.step(), LimitSide::right);
        eq_data_->cache_reallocate( *(eq_data_.get()), *(eq_data_.get()) );
        eq_data_->update_cache();
        auto p = *( eq_data_->bulk_int[0]->points(eq_data_->position_in_cache(eq_data_->computed_dh_cell_.elm_idx()), eq_data_.get()).begin() );

        arma::vec3 elm_cntr = this->eq_data_->computed_dh_cell_.elm().centre(); // {-0.5, 0.5, 0.0}

        for (uint i_comp=0; i_comp<3; ++i_comp) {
            EXPECT_EQ(elm_cntr(0), eq_data_->scalar_field[i_comp](p));
        }

        tg.next_time();
    }
}


TEST_F(MultiFieldTest, field_fe_test) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        vector_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: vector_fixed
          default_value: 0.0
    )YAML";
    std::vector< arma::vec3 > fe_expected = {{0.5, 0.5, -1.0}, {0.5, 1.0, 0.0}};

    this->read_input(eq_data_input);
    this->set_dh_cell(4, 39);

    for (uint i_time=0; i_time<2; i_time++) { // test in 2 time steps: 0.25, 1.0
        eq_data_->set_time(tg.step(), LimitSide::right);
        eq_data_->cache_reallocate( *(eq_data_.get()), *(eq_data_.get()) );
        eq_data_->update_cache();
        auto p = *( eq_data_->bulk_int[0]->points(eq_data_->position_in_cache(eq_data_->computed_dh_cell_.elm_idx()), eq_data_.get()).begin() );

        for (uint i_comp=0; i_comp<3; ++i_comp) {
            EXPECT_ARMA_EQ(fe_expected[i_time], eq_data_->vector_field[i_comp](p));
        }

        tg.next_time();
    }
}


TEST_F(MultiFieldTest, interpolated_p0_test) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          default_value: 0.0
          interpolation: P0_intersection
    )YAML";
    std::vector< double > p0_expected = {0.147606084199045, 1.147606084199045};

    this->read_input(eq_data_input);
    this->set_dh_cell(4, 39);

    for (uint i_time=0; i_time<2; i_time++) { // test in 2 time steps: 0.25, 1.0
        eq_data_->set_time(tg.step(), LimitSide::right);
        eq_data_->cache_reallocate( *(eq_data_.get()), *(eq_data_.get()) );
        eq_data_->update_cache();
        auto p = *( eq_data_->bulk_int[0]->points(eq_data_->position_in_cache(eq_data_->computed_dh_cell_.elm_idx()), eq_data_.get()).begin() );

        for (uint i_comp=0; i_comp<3; ++i_comp) {
            EXPECT_DOUBLE_EQ(p0_expected[i_time], eq_data_->scalar_field[i_comp](p));
        }

        tg.next_time();
    }
}



string eq_data_input = R"JSON(
[
  { id=37,
    a=1,
    b=0
  }
] 
}
)JSON";

TEST(Operators, assignment) {
    Profiler::instance();

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    auto mesh_reader = reader_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
    Mesh * mesh = mesh_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
    mesh_reader->read_physical_names(mesh);
    mesh_reader->read_raw_mesh(mesh);
    mesh->check_and_finish();

	std::vector<string> component_names = { "comp_0", "comp_1", "comp_2" };

	MultiField<3, FieldValue<3>::Scalar> mf_base;
	mf_base.name("a");
	mf_base.set_components(component_names);
	mf_base.units( UnitSI::dimensionless() );

	FieldSet set_of_field;
	set_of_field += mf_base;
    Input::Type::Array list_type = Input::Type::Array(set_of_field.make_field_descriptor_type("MultiFieldTest"));
    Input::ReaderToStorage reader( eq_data_input, list_type, Input::FileFormat::format_JSON);
    Input::Array in_list=reader.get_root_interface<Input::Array>();
    TimeGovernor tg(0.0, 1.0);
    set_of_field.set_input_list(in_list, tg);

    MultiField<3, FieldValue<3>::Scalar> mf_assignment;
	EXPECT_EQ("", mf_assignment.name());
	EXPECT_FALSE(mf_assignment.is_bc());

	// copies
	mf_assignment
	    .name("b")
	    .flags(FieldFlag::input_copy)
	    .units( UnitSI::dimensionless() );
	std::vector<string> component_names_2 = { "comp_a", "comp_b", "comp_c" };
	mf_assignment.set_components(component_names_2);
	mf_base.set_mesh( *mesh );
	mf_assignment.set_mesh( *mesh );
	mf_base.setup_components();

	MultiField<3, FieldValue<3>::Scalar> mf_copy(mf_base);	// copy constructor
	mf_assignment = mf_base; // assignment

	EXPECT_STREQ("a", mf_base.name().c_str());
	EXPECT_STREQ("b", mf_assignment.name().c_str());
	EXPECT_STREQ("a", mf_copy.name().c_str());
	EXPECT_EQ(3, mf_base.size());
	EXPECT_EQ(mf_assignment.size(), mf_base.size());
	EXPECT_EQ(mf_copy.size(), mf_base.size());
	for (unsigned int i=0; i<mf_base.size(); ++i) {
		EXPECT_EQ( component_names[i] + "_a", mf_base[i].name() );
		EXPECT_EQ( component_names_2[i] + "_b", mf_assignment[i].name() );
		EXPECT_EQ( mf_base[i].name(), mf_copy[i].name() );
	}

	{
		// throw assert, empty vector of component names
		MultiField<3, FieldValue<3>::Scalar> mf_assignment_error;
		mf_assignment_error
		    .name("c")
		    .flags(FieldFlag::input_copy);
		mf_assignment_error.set_mesh( *mesh );

		EXPECT_ASSERT_DEATH( { mf_assignment_error = mf_base; },
				"Vector of component names can't be empty!");
	}

	{
		// throw assert, source field has different component names
		std::vector<string> component_names = { "comp_a", "comp_b" };
		MultiField<3, FieldValue<3>::Scalar> mf_assignment_error;
		mf_assignment_error
		    .name("d")
		    .flags(FieldFlag::input_copy)
		    .units( UnitSI::dimensionless() );
		mf_assignment_error.set_components(component_names);
		mf_assignment_error.set_mesh( *mesh );
		mf_assignment_error.setup_components();

		EXPECT_ASSERT_DEATH( { mf_assignment_error = mf_base; },
				"Both multi fields must have same size of vectors");
	}

	delete mesh;
}
