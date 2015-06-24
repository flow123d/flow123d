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

#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>

#include <vector>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include "system/sys_profiler.hh"

#include "input/input_type.hh"
#include "input/type_output.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"

#include "fields/field_set.hh"
#include "fields/field_add_potential.hh"
#include "fields/unit_si.hh"
#include "fields/bc_field.hh"
#include "fields/multi_field.hh"
#include "coupling/equation.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/region.hh"
#include "flow/old_bcd.hh"
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

        template<int spacedim, class Value>
        class FieldFactory : public OldBcdInput::FieldFactory<spacedim, Value> {
        public:
        	FieldFactory( typename OldBcdInput::FieldFactory<spacedim, Value>::FieldPtr * field, bool read_flow = true )
        	: OldBcdInput::FieldFactory<spacedim, Value>(field),
        	  read_flow_(read_flow)
        	{}

        	virtual typename Field<spacedim,Value>::FieldBasePtr create_field(Input::Record rec, const FieldCommon &field) {
        		Input::AbstractRecord field_record;
        		if (rec.opt_val(field.input_name(), field_record)) {
        			return Field<spacedim,Value>::FieldBaseType::function_factory(field_record, field.n_comp() );
        		}
            	else {
            		OldBcdInput *old_bcd = OldBcdInput::instance();
            		if (read_flow_) {
            			old_bcd->read_flow_record(rec, field);
            		} else  {
            			old_bcd->read_transport_record(rec, field);
            		}
            		return *(this->field_);
            	}
        	}

        	/// Set true if field needs flow record, false if needs transport record
        	bool read_flow_;
        };

        EqData()
        {
        	arma::vec4 gravity = arma::vec4("3.0 2.0 1.0 -5.0");

            ADD_FIELD(anisotropy, "Anisothropic conductivity tensor.", "1.0");
            anisotropy.units( UnitSI::dimensionless() );

            ADD_FIELD(bc_type,"Boundary condition type, possible values:", "\"none\"" );
                      bc_type.input_selection(&get_bc_type_selection());
            bc_type.add_factory( std::make_shared<FieldFactory<3, FieldValue<3>::Enum> >
            					 (&(OldBcdInput::instance()->flow_type)) );
            bc_type.units( UnitSI::dimensionless() );

            ADD_FIELD(bc_pressure,"Dirichlet BC condition value for pressure." );
            bc_pressure.disable_where( bc_type, {none, neumann} );
        	bc_pressure.add_factory(std::make_shared< FieldAddPotential<3, FieldValue<3>::Scalar>::FieldFactory >
        			                (gravity, "bc_piezo_head"));

            bc_pressure.units( UnitSI::dimensionless() );

        	ADD_FIELD(bc_flux,"Flux in Neumman or Robin boundary condition." );
            bc_flux.disable_where( bc_type, {none, dirichlet, robin} );
        	bc_flux.add_factory( std::make_shared<FieldFactory<3, FieldValue<3>::Scalar> >
        						 (&(OldBcdInput::instance()->flow_flux)) );
            bc_flux.units( UnitSI::dimensionless() );

            ADD_FIELD(bc_robin_sigma,"Conductivity coefficient in Robin boundary condition.");
            bc_robin_sigma.disable_where( bc_type, {none, dirichlet, neumann} );
        	bc_robin_sigma.add_factory( std::make_shared<FieldFactory<3, FieldValue<3>::Scalar> >
        								(&(OldBcdInput::instance()->flow_sigma)) );
            bc_robin_sigma.units( UnitSI::dimensionless() );

            ADD_FIELD(bc_conc, "BC concentration", "0.0" );
            bc_conc.add_factory( std::make_shared<FieldFactory<3, FieldValue<3>::Vector> >
            					 (&(OldBcdInput::instance()->trans_conc), false) );
            bc_conc.units( UnitSI::dimensionless() );
        }

        Field<3, FieldValue<3>::TensorFixed > anisotropy;
        BCField<3, FieldValue<3>::Enum > bc_type; // Discrete need Selection for initialization
        BCField<3, FieldValue<3>::Scalar > bc_pressure;
        BCField<3, FieldValue<3>::Scalar > bc_flux;
        BCField<3, FieldValue<3>::Scalar > bc_robin_sigma;
        BCField<3, FieldValue<3>::Vector > bc_conc;

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
        Field<3, FieldValue<3>::Vector > init_conc;
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
        component_names = { "comp_0", "comp_1", "comp_2", "comp_3" };

    }

    void read_input(const string &input) {
        // read input string
        Input::JSONToStorage reader( input, get_input_type() );
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

        DBGMSG("init\n");
        data.set_mesh(*mesh);
        data.set_input_list( in_rec.val<Input::Array>("data"));
        data.set_limit_side(LimitSide::right);
        data.set_time(tg.step());
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
	        		IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_decsription("SomeEquation_Data") )
	                .copy_keys( SomeEquation::EqData().make_field_descriptor_type("SomeEquation") )
	                .declare_key("bc_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type(nullptr), "" )
	                .declare_key(OldBcdInput::flow_old_bcd_file_key(), IT::FileName::input(), "")
	                .declare_key(OldBcdInput::transport_old_bcd_file_key(), IT::FileName::input(), "")
	                .declare_key("init_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type(nullptr), "" )
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
            init_conc = [ 1, 2, 3, 4],
            conc_mobile = {
                TYPE="MultiField",
                component_names=["comp_0", "comp_1", "comp_2", "comp_3"],
                common={TYPE="FieldConstant", value=[1, 2, 3, 4]},
                components=[ {TYPE="FieldConstant", value=1}, {TYPE="FieldConstant", value=2}, {TYPE="FieldConstant", value=3}, {TYPE="FieldConstant", value=4}]
              }
          },
          { region="2D XY diagonal",
            init_pressure=2.2,
            conc_mobile={REF="/data/0/conc_mobile"}
          },
          { r_set="BULK",
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
            conc_mobile = {
                TYPE="MultiField",
                component_names=["comp_0", "comp_1", "comp_2", "comp_3"],
                common={TYPE="FieldConstant", value=[5, 6, 7, 8]},
                components=[ {TYPE="FieldConstant", value=5}, {TYPE="FieldConstant", value=6}, {TYPE="FieldConstant", value=7}, {TYPE="FieldConstant", value=8}]
              }
          },
          { rid=102,
            bc_type="dirichlet",
            bc_piezo_head=1.23,
            conc_mobile={REF="/data/3/conc_mobile"}
          }
      ] 
    }
    )JSON";

    read_input(eq_data_input);
    // cout << Input::Type::OutputText(&SomeEquation::get_input_type()) << endl;

    Space<3>::Point p;
    p(0)=1.0; p(1)= 2.0; p(2)=3.0;

    DBGMSG("elements size: %d %d\n",mesh->element.size(), mesh->bc_elements.size());

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
    for (unsigned int i=0; i<data.conc_mobile.size(); ++i) {     // multifield
        EXPECT_DOUBLE_EQ( 1.0 + i, data.conc_mobile[i].value(p, el_1d) );
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
    for (unsigned int i=0; i<data.conc_mobile.size(); ++i) {     // multifield
        EXPECT_DOUBLE_EQ( 1.0 + i, data.conc_mobile[i].value(p, el_2d) );
    }

    // init_conc - variable length vector
    FieldValue<3>::Vector::return_type conc = data.init_conc.value(p, el_1d);
    EXPECT_EQ(1 ,conc.n_cols);
    EXPECT_EQ(4 ,conc.n_rows);
    EXPECT_DOUBLE_EQ(1 ,conc[0]);
    EXPECT_DOUBLE_EQ(2 ,conc[1]);
    EXPECT_DOUBLE_EQ(3 ,conc[2]);
    EXPECT_DOUBLE_EQ(4 ,conc[3]);

    // bulk_set_filed - test setting on region set, test setting field without default value
    EXPECT_EQ( 5.7, data.bulk_set_field.value(p, el_1d) );
    EXPECT_EQ( 5.7, data.bulk_set_field.value(p, el_2d) );
    EXPECT_EQ( 5.7, data.bulk_set_field.value(p, el_3d) );

    //boundary fields
    EXPECT_EQ( EqData::dirichlet, data.bc_type.value(p, el_bc_top) );
    EXPECT_DOUBLE_EQ(1.23, data.bc_pressure.value(p, el_bc_top) );    // pressure

    EXPECT_EQ( EqData::dirichlet, data.bc_type.value(p, el_bc_bottom) );
    EXPECT_DOUBLE_EQ(1.23 + (3 + 4 + 3 - 5), data.bc_pressure.value(p, el_bc_bottom) );    // piezo_head

}



TEST_F(SomeEquation, old_bcd_input) {
    // Test input for old_bcd
    string eq_data_old_bcd = R"JSON(
    { 
      data=[
          { r_set="BOUNDARY",
            flow_old_bcd_file="coupling/simplest_cube.fbc",
            transport_old_bcd_file="coupling/transport.fbc",
            conc_mobile = {
                TYPE="MultiField",
                component_names=["comp_0", "comp_1", "comp_2", "comp_3"],
                common={TYPE="FieldConstant", value=[1, 2, 3, 4]},
                components=[ {TYPE="FieldConstant", value=1}, {TYPE="FieldConstant", value=2}, {TYPE="FieldConstant", value=3}, {TYPE="FieldConstant", value=4}]
              }
          },
          { r_set="BULK",
            bulk_set_field=0.0,
            conc_mobile={REF="/data/0/conc_mobile"}
          } 
      ] 
    }
    )JSON";

    read_input(eq_data_old_bcd);

    Space<3>::Point p;
    p(0)=1.0; p(1)= 2.0; p(2)=3.0;

    DBGMSG("elements size: %d %d\n",mesh->element.size(), mesh->bc_elements.size());


    // Four bc elements are read with mesh, corresponding to BCD IDs:
    // 7, 17, 10, 12
    // The bcd IDs  order in the bc_vector: 7, 17, 10, 12, 0, 1, 2, 3, 4, 5, 6, 8, 9, 11, 13, 14, 15, 16

    {
    auto test_elm = mesh->element_accessor(4, true);
    EXPECT_EQ(EqData::dirichlet, (EqData::BC_type)data.bc_type.value(p, test_elm));
    EXPECT_DOUBLE_EQ(1.0, data.bc_pressure.value(p, test_elm) );
    arma::vec value = data.bc_conc.value(p, test_elm);
    EXPECT_DOUBLE_EQ(1.0, value(0) );
    EXPECT_DOUBLE_EQ(11.0, value(1) );
    EXPECT_DOUBLE_EQ(21.0, value(2) );
    EXPECT_DOUBLE_EQ(31.0, value(3) );
    }

    {
    auto test_elm = mesh->element_accessor(10, true);
    EXPECT_EQ(EqData::dirichlet, (EqData::BC_type)data.bc_type.value(p, test_elm));
    EXPECT_DOUBLE_EQ(7.0, data.bc_pressure.value(p, test_elm) );
    arma::vec value = data.bc_conc.value(p, test_elm);
    EXPECT_DOUBLE_EQ(7.0, value(0) );
    EXPECT_DOUBLE_EQ(17.0, value(1) );
    EXPECT_DOUBLE_EQ(27.0, value(2) );
    EXPECT_DOUBLE_EQ(37.0, value(3) );
    }

}



TEST_F(SomeEquation, wrong_time_order) {
    // Test input for old_bcd
    string eq_data = R"JSON(
    { 
      data=[
          { r_set="BULK",
            time=1.0,
            bulk_set_field=0.0
          }, 
          { r_set="BULK",
            time=0.0,
            bulk_set_field=1.0
          } 
      ] 
    }
    )JSON";

    EXPECT_THROW({read_input(eq_data);}, FieldCommon::ExcNonascendingTime);
}
