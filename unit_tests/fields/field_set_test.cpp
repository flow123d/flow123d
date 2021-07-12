/*
 * filed_set_test_.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include "fields/field_set.hh"
#include "tools/unit_si.hh"
#include "fields/bc_field.hh"
#include "fields/field_model.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"


FLOW123D_FORCE_LINK_IN_PARENT(field_constant)
FLOW123D_FORCE_LINK_IN_PARENT(field_formula)


enum {
	r_first,
	r_second,
	r_third
};

auto reaction_type_sel = Input::Type::Selection("ReactionType")
		.add_value(r_first, "r_first")
		.add_value(r_second, "r_second")
		.add_value(r_third, "r_third")
		.close();

// Test input for 'values' test
const string eq_data_input = R"JSON(
[
      { 
        time=0.0,
        region="BULK",
        init_pressure=1.1,
        velocity={TYPE="FieldFormula",
            value=[ "x", "y" , "z"]
        },
        reaction_type="r_first"
      },
      { 
        time=1.0,
        region="BULK",
        velocity=[1,2,4]
      }
]
)JSON";

class SomeEquation : public testing::Test {

public:
	class EqData : public FieldSet {
	public:
		EqData() {
			*this += velocity
						.name("velocity")
						.description("Velocity vector.")
						.input_default("0.0")
						.flags_add(in_main_matrix)
						.units( UnitSI().kg(3).m() );
			*this += init_pressure
						.disable_where(type, {r_first, r_second })
						.name("init_pressure")
						.description("Pressure head")
						.units( UnitSI().m() );

			*this += type
						.name("reaction_type")
						.description("")
						.units( UnitSI::dimensionless() )
						.flags_add(in_main_matrix)
						.input_selection(reaction_type_sel);
		}

		// fields
	    Field<3, FieldValue<3>::VectorFixed > velocity;
	    Field<3, FieldValue<3>::Scalar > init_pressure;
	    Field<3, FieldValue<3>::Enum > type;
	};

	SomeEquation() {
	    Profiler::instance();
	    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        mesh_ = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
	}

	~SomeEquation() {
		delete mesh_;
	}

	Mesh * mesh_;
	std::vector<string> component_names_;
};


TEST_F(SomeEquation, add_operator_death) {
    auto data = EqData();
    Field<3, FieldValue<3>::Scalar > pressure;
    EXPECT_ASSERT_DEATH({
        data+=pressure
             .name("init_pressure");},
          "Another field of the same name exists");

}


TEST_F(SomeEquation, access) {
	auto data = EqData();

	// test basic add in EqData constructor
	// test size method
	// test access operator []
	EXPECT_EQ(3, data.size() );

	EXPECT_EQ(&data.velocity,  &(data["velocity"]) );
    EXPECT_EQ("velocity",  data["velocity"].name() );
    EXPECT_EQ("Velocity vector.",  data["velocity"].description() );
	EXPECT_EQ(&data.init_pressure,  &(data["init_pressure"]) );
	EXPECT_EQ(&data.type,  &(data["reaction_type"]) );

	// test subset
	FieldSet set=data.subset({"velocity", "init_pressure"});

	EXPECT_EQ(2, set.size());
	EXPECT_EQ(&data.velocity,  &(set["velocity"]) );
    EXPECT_EQ(&data.init_pressure,  &(set["init_pressure"]) );
    EXPECT_EQ(nullptr, set.field("reaction_type"));

    // test add of other field set
    Field<3, FieldValue<3>::Scalar > pressure;
    Field<3, FieldValue<3>::Scalar > copy_pressure;
    FieldSet other;

    other+=pressure.name("pressure");
    other+=copy_pressure.name("copy_pressure");
    data+=other;
    EXPECT_EQ(5, data.size());
    EXPECT_EQ(&pressure,  &(data["pressure"]) );
    EXPECT_EQ(&copy_pressure,  &(data["copy_pressure"]) );
}

TEST_F(SomeEquation, access_death) {
    auto data = EqData();

    EXPECT_THROW({data["noname"];}, FieldSet::ExcUnknownField);

}

TEST_F(SomeEquation, subset_mask) {
    auto data = EqData();

    auto main_matrix_set=data.subset(data.in_main_matrix);
    EXPECT_EQ(2, main_matrix_set.size());
    EXPECT_EQ(&data.velocity,  &(main_matrix_set["velocity"]) );
    EXPECT_EQ(&data.type,  &(main_matrix_set["reaction_type"]) );
    EXPECT_EQ(nullptr, main_matrix_set.field("init_pressure") );
}


TEST_F(SomeEquation, field_descriptor) {
	Input::Type::Record descriptor = EqData().make_field_descriptor_type("SomeEquation");

	descriptor.finish();
	EXPECT_EQ(7, descriptor.size());
	EXPECT_TRUE( descriptor.has_key("time"));
	EXPECT_TRUE( descriptor.has_key("rid"));
	EXPECT_TRUE( descriptor.has_key("region"));
	EXPECT_TRUE( descriptor.has_key("velocity"));
	EXPECT_TRUE( descriptor.has_key("init_pressure"));
	EXPECT_TRUE( descriptor.has_key("reaction_type"));
}



TEST_F(SomeEquation, set_field) {
    auto data = EqData();
    Field<3, FieldValue<3>::Scalar > other_pressure;
    other_pressure
    .name("other_pressure")
    .description("other pressure");
    other_pressure.set_mesh(*mesh_);

    EXPECT_EQ("init_pressure", data["init_pressure"].name());
    EXPECT_EQ("Pressure head", data["init_pressure"].description());
    data["init_pressure"].flags(FieldFlag::input_copy);
    data.set_field("init_pressure", other_pressure);
    EXPECT_EQ(nullptr, data.field("other_pressure"));
    EXPECT_EQ("init_pressure", data["init_pressure"].name());
    EXPECT_EQ("other pressure", data["init_pressure"].description());

}


TEST_F(SomeEquation, collective_interface) {
    auto data = EqData();
    component_names_ = { "component_0", "component_1", "component_2", "component_3" };

    data.set_components(component_names_);
    EXPECT_EQ(component_names_.size(), data["init_pressure"].n_comp());
    EXPECT_EQ(component_names_.size(), data["velocity"].n_comp());
    EXPECT_EQ(component_names_.size(), data["reaction_type"].n_comp());

    EXPECT_EQ(nullptr,data["init_pressure"].mesh());
    data.set_mesh(*mesh_);
    EXPECT_EQ(mesh_, data["init_pressure"].mesh());
    EXPECT_EQ(mesh_, data["velocity"].mesh());
    EXPECT_EQ(mesh_, data["reaction_type"].mesh());

    // flags_add
    FieldFlag::Flags matrix(
        FieldSet::equation_input
        & FieldSet::declare_input
        & FieldSet::allow_output
        & FieldSet::in_main_matrix
        & FieldSet::in_rhs);
    FieldFlag::Flags rhs(
    FieldSet::equation_input
      & FieldSet::declare_input
      & FieldSet::allow_output
      & FieldSet::in_rhs);

    data.flags_add(FieldSet::in_rhs);
    EXPECT_EQ( rhs, data["init_pressure"].flags() );
    EXPECT_EQ( matrix, data["velocity"].flags() );
    EXPECT_EQ( matrix, data["reaction_type"].flags() );

    auto output_types = OutputTimeSet::empty_discrete_flags();
    output_types[OutputTime::CORNER_DATA] = true;
    data.output_type(output_types);
    EXPECT_TRUE( data["init_pressure"].get_output_type()[OutputTime::CORNER_DATA] );
    EXPECT_TRUE( data["velocity"].get_output_type()[OutputTime::CORNER_DATA] );
    EXPECT_TRUE( data["reaction_type"].get_output_type()[OutputTime::CORNER_DATA] );
}

TEST_F(SomeEquation, input_related) {
    auto data = EqData();
    component_names_ = { "component_0", "component_1" };

    data.set_components(component_names_);
    Input::Type::Array list_type = Input::Type::Array(data.make_field_descriptor_type("SomeEquation"));
    Input::ReaderToStorage reader( eq_data_input, list_type, Input::FileFormat::format_JSON);
    Input::Array in=reader.get_root_interface<Input::Array>();
    TimeGovernor tg(0.0, 0.5);
    data.set_input_list(in, tg);
    data.set_mesh(*mesh_);

    data.mark_input_times(tg);
    Region front_3d = mesh_->region_db().find_id(40);
    // time = 0.0
    data.set_time(tg.step(), LimitSide::right);
    EXPECT_FALSE(data.is_constant(front_3d));
    EXPECT_TRUE(data.changed());
    EXPECT_TRUE(tg.is_current(tg.marks().type_input()));
    data.set_time(tg.step(), LimitSide::right);
    EXPECT_TRUE(data.changed());
    tg.next_time();

    // time = 0.5
    data.set_time(tg.step(), LimitSide::right);
    EXPECT_FALSE(data.changed());
    EXPECT_FALSE(data.is_constant(front_3d));
    EXPECT_FALSE(tg.is_current(tg.marks().type_input()));
    tg.next_time();

    // time = 1.0
    data.set_time(tg.step(), LimitSide::right);
    EXPECT_TRUE(data.changed());
    EXPECT_TRUE(data.is_constant(front_3d));
    EXPECT_TRUE(tg.is_current(tg.marks().type_input()));

}


/******************************
 *  TEST OF FIELD DEPENDENCY  *
 ******************************/
namespace IT = Input::Type;

// Functor computing FieldModels
struct fn_f_field {
    inline arma::vec3 operator() (double a, arma::vec3 v) {
        return a * v;
    }
};
struct fn_c_field {
    inline double operator() (arma::vec3 v) {
        return arma::norm(v, 2);
    }
};
struct fn_d_field {
    inline double operator() (double a, double b) {
        return a * b;
    }
};
struct fn_g_field {
    inline arma::vec3 operator() (double a, arma::vec3 u, arma::vec3 v) {
        return u + a * v;
    }
};
struct fn_e_field {
    inline double operator() (double a) {
        return 0.5 * a;
    }
};

class TestDependency : public testing::Test {

public:
    class EqData : public FieldSet {
    public:
        EqData() {
            *this += e_field
                    .name("e_field")
                    .description("Model: 0.5 * d_field")
                    .units( UnitSI().m() );
            *this += d_field
                    .name("d_field")
                    .description("Model: a_field * c_field")
                    .units( UnitSI::dimensionless() );
            *this += c_field
                    .name("c_field")
                    .description("Model: norm(b_field)")
                    .units( UnitSI().m() );
            *this += b_field
                    .name("b_field")
                    .description("Vector field")
                    .units( UnitSI().m() );
            *this += g_field
                    .name("g_field")
                    .description("Model: f_field + d_field * b_field")
                    .units( UnitSI::dimensionless() );
            *this += f_field
                    .name("f_field")
                    .description("Model: a_field * b_field")
                    .units( UnitSI::dimensionless() );
            *this += a_field
                    .name("a_field")
                    .description("Scalar field.")
                    .units( UnitSI::dimensionless() );
        }

        /// Initialize FieldModels
        void initialize() {
        	c_field.set(Model<3, FieldValue<3>::Scalar>::create(fn_c_field(), b_field), 0.0);
            d_field.set(Model<3, FieldValue<3>::Scalar>::create(fn_d_field(), a_field, c_field), 0.0);
            e_field.set(Model<3, FieldValue<3>::Scalar>::create(fn_e_field(), d_field), 0.0);
            f_field.set(Model<3, FieldValue<3>::VectorFixed>::create(fn_f_field(), a_field, b_field), 0.0);
            g_field.set(Model<3, FieldValue<3>::VectorFixed>::create(fn_g_field(), d_field, f_field, b_field), 0.0);
        }

        /// Check result
        void check() {
        	/* print
            std::cout << "Original order in field_list:\n";
            for (auto f : this->field_list) std::cout << " " << f->name();
            std::cout << "\n----------------------------------------------------\n";
            std::cout << "Sorted order by regions:\n";
            for (auto r : this->region_field_update_order_) {
                std::cout << " " << r.first << ":  ";
                for (auto f : r.second) std::cout << " " << f->name();
                std::cout << "\n";
            }
            std::cout << "----------------------------------------------------\n";
            // */

            // check
            std::vector<std::string> orig_order; // holds original order in field_list
            for (auto f : this->field_list) orig_order.push_back(f->name());
            for (auto r : this->region_field_update_order_) {
                EXPECT_EQ(r.second.size(), 7);
                if (r.first % 2) // only bulk regions are sorted
                    for (unsigned int i=1; i<r.second.size(); ++i)
                        // fields are sorted by name
                        EXPECT_TRUE(r.second[i-1]->name() < r.second[i]->name());
                else // boundary regions are in original order
                	for (unsigned int i=0; i<r.second.size(); ++i)
                	    EXPECT_EQ(r.second[i]->name(), orig_order[i]);
            }
        }

        // fields
        Field<3, FieldValue<3>::Scalar >      a_field;
        Field<3, FieldValue<3>::VectorFixed > b_field;
        Field<3, FieldValue<3>::Scalar >      c_field;
        Field<3, FieldValue<3>::Scalar >      d_field;
        Field<3, FieldValue<3>::Scalar >      e_field;
        Field<3, FieldValue<3>::VectorFixed > f_field;
        Field<3, FieldValue<3>::VectorFixed > g_field;
    };

    TestDependency() {
        Profiler::instance();
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        mesh_ = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
        data_ = std::make_shared<EqData>();
    }

    ~TestDependency() {
        delete mesh_;
    }

    void read_input(const string &input) {
        // read input string
        Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_YAML );
        Input::Record in_rec=reader.get_root_interface<Input::Record>();

        TimeGovernor tg(0.0, 1.0);

        //data.set_components(component_names);        // set number of substances posibly read from elsewhere

        static std::vector<Input::Array> inputs;
        unsigned int input_last = inputs.size(); // position of new item
        inputs.push_back( in_rec.val<Input::Array>("data") );

        data_->set_mesh(*mesh_);
        data_->initialize();
        data_->set_input_list( inputs[input_last], tg );
        data_->set_time(tg.step(), LimitSide::right);
        data_->set_dependency( *(data_.get()) );
    }

    static Input::Type::Record & get_input_type() {
        return IT::Record("TestDependency","")
                .declare_key("data", IT::Array(
                        IT::Record("TestDependency_Data", FieldCommon::field_descriptor_record_description("TestDependency_Data") )
                        .copy_keys( TestDependency::EqData().make_field_descriptor_type("TestDependency") )
                        .close()
                        ), IT::Default::obligatory(), ""  )
                .close();
    }

    Mesh * mesh_;
    std::shared_ptr<EqData> data_;
};

string dependency_input = R"YAML(
data:
  - region: ALL
    time: 0.0
    a_field: !FieldConstant
      value: 0.5
    b_field: [1, 2, 3]
)YAML";



TEST_F(TestDependency, dependency_tree) {
	this->read_input(dependency_input);
    data_->check();
}


    /*
     * set_time
     */


    /*
     * output
     */



