/*
 * eq_data_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
 *
 *
 *  This unit test is meant as proof of concept for EqDataBase + Fields.
 *  We demonstrate and test:
 *  - inheritance of fields from some base class that do not perform any action with the data
 *  - declaration of own fields and their registration in EqData constructor
 *  - reading fields from the input
 *  - accessing field data
 */

/*
 * TODO:
 * - Meli bychom mit vsechny rovnice pojmenovane? Jake schema mit pro pojmenovani recordu EqData, ma mit kazda instance konkretni rovnice vlastni
 *   Record (nesmysl)? Tedy v nazvu recordu pouzit nazev tridy
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include <vector>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include "system/sys_profiler.hh"

#include "input/input_type.hh"
#include "input/type_output.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"

#include "fields/field_set.hh"
#include "fields/field_add_potential.hh"
#include "fields/unit_si.hh"
#include "fields/bc_field.hh"
#include "fields/multi_field.hh"
#include "coupling/equation.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/region.hh"
#include <armadillo>

using namespace std;
namespace IT=Input::Type;





class SomeEquationBase : public EquationBase {
protected:
    class EqData : public FieldSet {
    public:
        enum BC_type {
            none=0,
            dirichlet=1,
            neumann=2,
            robin=3,
            total_flux=4
        };

        static const IT::Selection & get_bc_type_selection();

        EqData()
        {
        	arma::vec4 gravity = arma::vec4("3.0 2.0 1.0 -5.0");

            ADD_FIELD(anisotropy, "Anisothropic conductivity tensor.", "1.0");
            anisotropy.units( UnitSI::dimensionless() );

            ADD_FIELD(bc_type,"Boundary condition type, possible values:", "\"none\"" );
                      bc_type.input_selection(get_bc_type_selection());
            bc_type.units( UnitSI::dimensionless() );

            ADD_FIELD(bc_pressure,"Dirichlet BC condition value for pressure." );
            bc_pressure.disable_where( bc_type, {none, neumann} );
        	bc_pressure.add_factory(std::make_shared< FieldAddPotential<3, FieldValue<3>::Scalar>::FieldFactory >
        			                (gravity, "bc_piezo_head"));

            bc_pressure.units( UnitSI::dimensionless() );

        	ADD_FIELD(bc_flux,"Flux in Neumman or Robin boundary condition." );
            bc_flux.disable_where( bc_type, {none, dirichlet, robin} );
            bc_flux.units( UnitSI::dimensionless() );

            ADD_FIELD(bc_robin_sigma,"Conductivity coefficient in Robin boundary condition.");
            bc_robin_sigma.disable_where( bc_type, {none, dirichlet, neumann} );
            bc_robin_sigma.units( UnitSI::dimensionless() );

            ADD_FIELD(bc_conc, "BC concentration", "0.0" );
            bc_conc.units( UnitSI::dimensionless() );
        }

        Field<3, FieldValue<3>::TensorFixed > anisotropy;
        BCField<3, FieldValue<3>::Enum > bc_type; // Discrete need Selection for initialization
        BCField<3, FieldValue<3>::Scalar > bc_pressure;
        BCField<3, FieldValue<3>::Scalar > bc_flux;
        BCField<3, FieldValue<3>::Scalar > bc_robin_sigma;
        BCField<3, FieldValue<3>::VectorFixed > bc_conc;

    };

    void output_data() override {}
};

const IT::Selection & SomeEquationBase::EqData::get_bc_type_selection() {
	return IT::Selection("EqData_bc_Type")
             .add_value(none, "none")
             .add_value(dirichlet, "dirichlet")
             .add_value(neumann, "neumann")
             .add_value(robin, "robin")
             .add_value(total_flux, "total_flux")
			 .close();
}



class SomeEquation : public testing::Test, SomeEquationBase {
public:

    class EqData : public SomeEquationBase::EqData {
    public:

        EqData() : SomeEquationBase::EqData() {
            ADD_FIELD(init_pressure, "Initial condition as pressure", "0.0" );
            ADD_FIELD(init_conc, "Initial condition for the concentration (vector of size equal to n. components", "0.0" );
            ADD_FIELD(bulk_set_field, "");
            ADD_FIELD(conc_mobile, "");

            init_pressure.units( UnitSI::dimensionless() );
            init_conc.units( UnitSI::dimensionless() );
            bulk_set_field.units( UnitSI::dimensionless() );
            conc_mobile.units( UnitSI::dimensionless() );
        }

        Field<3, FieldValue<3>::Scalar > init_pressure;
        Field<3, FieldValue<3>::VectorFixed > init_conc;
        Field<3, FieldValue<3>::Scalar > bulk_set_field;
        MultiField<3, FieldValue<3>::Scalar > conc_mobile;
    };

protected:
    static const Input::Type::Record & get_input_type();
    static MultiField<3, FieldValue<3>::Scalar> empty_mf;
    EqData data;
    std::vector<string> component_names;

    virtual void SetUp() {
        Profiler::initialize();
        //data.gravity_=arma::vec4("3.0 2.0 1.0 -5.0");
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

        FilePath mesh_file("mesh/simplest_cube.msh", FilePath::input_file);
        mesh= new Mesh;
        ifstream in(string( mesh_file ).c_str());
        mesh->read_gmsh_from_stream(in);
        component_names = { "comp_0", "comp_1", "comp_2" };

    }

    void read_input(const string &input) {
        // read input string
        Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_JSON );
        Input::Record in_rec=reader.get_root_interface<Input::Record>();

        TimeGovernor tg(0.0, 1.0);

        data.set_components(component_names);        // set number of substances posibly read from elsewhere
        //data.bc_conc.set_components(component_names);

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

        DebugOut() << "init\n";

        static std::vector<Input::Array> inputs;
        unsigned int input_last = inputs.size(); // position of new item
        inputs.push_back( in_rec.val<Input::Array>("data") );

        data.set_mesh(*mesh);
        data.set_input_list( inputs[input_last] );
        data.set_time(tg.step(), LimitSide::right);
    }

    virtual void TearDown() {
        delete mesh;
    };


    Mesh *mesh;
};



MultiField<3, FieldValue<3>::Scalar> SomeEquation::empty_mf = MultiField<3, FieldValue<3>::Scalar>();

const IT::Record & SomeEquation::get_input_type() {
	return IT::Record("SomeEquation","")
	        .declare_key("data", IT::Array(
	        		IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
	                .copy_keys( SomeEquation::EqData().make_field_descriptor_type("SomeEquation") )
	                .declare_key("bc_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
	                .declare_key("init_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
					.close()
	                ), IT::Default::obligatory(), ""  )
			.close();
}



TEST_F(SomeEquation, values) {
    // Test input for 'values' test
    string eq_data_input = R"JSON(
    { 
      data=[
          { rid=37,
            init_pressure={
                TYPE="FieldConstant",
                value=1.1
              },
            init_conc = [ 1, 2, 3 ],
            // MultiField
            conc_mobile = {
                TYPE="FieldConstant", 
                value=[1, 2, 3]
              }
          },
          { region= ["2D XY diagonal", "3D back"],
            init_pressure=2.2,
            conc_mobile={REF="/data/0/conc_mobile"}
          },
          { region="BULK",
            bulk_set_field=5.7,
            conc_mobile={REF="/data/0/conc_mobile"}
          },

          // boundary          
          { rid=101,
            bc_type={TYPE="FieldConstant", value = "dirichlet"},
            bc_pressure={
                TYPE="FieldConstant",
                value=1.23
            },
            bc_conc={
                TYPE="FieldFormula",
                value=["x", "10+x", "20+x"]
            },
            conc_mobile = {
                TYPE="FieldConstant", 
                value=[5, 6, 7]
              }
          },
          { rid=102,
            bc_type="dirichlet",
            bc_piezo_head=1.23,
            bc_conc={REF="/data/3/bc_conc"},
            conc_mobile={REF="/data/3/conc_mobile"}
          }
      ] 
    }
    )JSON";

    read_input(eq_data_input);

    Space<3>::Point p;
    p(0)=1.0; p(1)= 2.0; p(2)=3.0;

    DebugOut().fmt("elements size: {} {}\n", mesh->element.size(), mesh->bc_elements.size());

    // check element accessors
    ElementAccessor<3> el_1d=mesh->element_accessor(0); // region 37 "1D diagonal"
    EXPECT_EQ(37, el_1d.region().id());
    ElementAccessor<3> el_2d=mesh->element_accessor(1); // region 38 "2D XY diagonal"
    EXPECT_EQ(38, el_2d.region().id());
    ElementAccessor<3> el_3d=mesh->element_accessor(3); // region 39 "3D back"
    EXPECT_EQ(39, el_3d.region().id());
    ElementAccessor<3> el_bc_top=mesh->element_accessor(0,true); // region 101 ".top side"
    EXPECT_EQ(101, el_bc_top.region().id());
    ElementAccessor<3> el_bc_bottom=mesh->element_accessor(2,true); // region 102 ".top side"
    EXPECT_EQ(102, el_bc_bottom.region().id());

    // bulk fields
    EXPECT_DOUBLE_EQ(1.1, data.init_pressure.value(p, el_1d) );
    auto conc_mobile_val = data.conc_mobile.value(p, el_1d);
    for (unsigned int i=0; i<data.conc_mobile.size(); ++i) {     // multifield
        EXPECT_DOUBLE_EQ( 1.0 + i, conc_mobile_val[i] );
    }

    FieldValue<3>::TensorFixed::return_type value = data.anisotropy.value(p, el_1d);
    EXPECT_DOUBLE_EQ( 1.0, value.at(0,0) );
    EXPECT_DOUBLE_EQ( 0.0, value.at(0,1) );
    EXPECT_DOUBLE_EQ( 0.0, value.at(0,2) );

    EXPECT_DOUBLE_EQ( 0.0, value.at(1,0) );
    EXPECT_DOUBLE_EQ( 1.0, value.at(1,1) );
    EXPECT_DOUBLE_EQ( 0.0, value.at(1,2) );

    EXPECT_DOUBLE_EQ( 0.0, value.at(2,0) );
    EXPECT_DOUBLE_EQ( 0.0, value.at(2,1) );
    EXPECT_DOUBLE_EQ( 1.0, value.at(2,2) );

    EXPECT_DOUBLE_EQ(2.2, data.init_pressure.value(p, el_2d) );
    EXPECT_DOUBLE_EQ(2.2, data.init_pressure.value(p, el_3d) );
    conc_mobile_val = data.conc_mobile.value(p, el_2d);
    for (unsigned int i=0; i<data.conc_mobile.size(); ++i) {     // multifield
        EXPECT_DOUBLE_EQ( 1.0 + i, conc_mobile_val[i] );
    }

    // init_conc - fixed length vector
    FieldValue<3>::VectorFixed::return_type conc = data.init_conc.value(p, el_1d);
    EXPECT_DOUBLE_EQ(1 ,conc[0]);
    EXPECT_DOUBLE_EQ(2 ,conc[1]);
    EXPECT_DOUBLE_EQ(3 ,conc[2]);

    // bulk_set_filed - test setting on region set, test setting field without default value
    EXPECT_EQ( 5.7, data.bulk_set_field.value(p, el_1d) );
    EXPECT_EQ( 5.7, data.bulk_set_field.value(p, el_2d) );
    EXPECT_EQ( 5.7, data.bulk_set_field.value(p, el_3d) );

    //boundary fields
    EXPECT_EQ( EqData::dirichlet, data.bc_type.value(p, el_bc_top) );
    EXPECT_DOUBLE_EQ(1.23, data.bc_pressure.value(p, el_bc_top) );    // pressure

    EXPECT_EQ( EqData::dirichlet, data.bc_type.value(p, el_bc_bottom) );
    EXPECT_DOUBLE_EQ(1.23 + (3 + 4 + 3 - 5), data.bc_pressure.value(p, el_bc_bottom) );    // piezo_head

    arma::vec bc_value = data.bc_conc.value(p, el_bc_bottom);
    EXPECT_DOUBLE_EQ(1.0, bc_value(0) );
    EXPECT_DOUBLE_EQ(11.0, bc_value(1) );
    EXPECT_DOUBLE_EQ(21.0, bc_value(2) );
}



TEST_F(SomeEquation, wrong_time_order) {
    // Test input for wrong_time_order
    string eq_data = R"JSON(
    { 
      data=[
          { region="BULK",
            time=1.0,
            bulk_set_field=0.0
          }, 
          { region="BULK",
            time=0.0,
            bulk_set_field=1.0
          } 
      ] 
    }
    )JSON";

    EXPECT_THROW({read_input(eq_data);}, FieldCommon::ExcNonascendingTime);
}
