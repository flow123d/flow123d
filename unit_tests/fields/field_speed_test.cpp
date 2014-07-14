/*
 * field_speed_test.cpp
 *
 *  Created on: Feb 3, 2013
 *      Author: jb
 */

#define TEST_USE_MPI

#include <flow_gtest_mpi.hh>

#include "fields/field_constant.hh"
#include "fields/field.hh"
#include "fields/field_set.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"


string input = R"INPUT(
{   
    scalar_value=[1.1, 2.2, 3.3],
    vector_fixed_value=[1.2, 2.3, 3.4],
    vector_value=[4.5, 5.6, 6.7],
    tensor_fixed_value=[ [1.2, 2.3], [3.4, 4.5]],

    scalar_data=[
	  { rid=37,
		init_pressure={
			TYPE="FieldConstant",
			value=1.1
		  }
	  },
	  { rid=101,
		init_pressure={
			TYPE="FieldConstant",
			value=2.2
		  }
	  },
	  { rid=102,
		init_pressure={
			TYPE="FieldConstant",
			value=3.3
		  }
	  }
    ],
    vector_fixed={
        TYPE="FieldConstant",
        value=[1.2, 2.3, 3.4]
    },
    vector={
        TYPE="FieldConstant",
        value=[4.5, 5.6, 6.7]
    },
    tensor_fixed={
        TYPE="FieldConstant",
        value=[ 1.2, 2.3, 3.4, 4.5 ]
    }
}
)INPUT";


class FieldSpeed : public testing::Test {
public:
	typedef FieldConstant<3, FieldValue<3>::Scalar > ScalarField;
    typedef FieldConstant<3, FieldValue<3>::VectorFixed > VecFixField;
    typedef FieldConstant<3, FieldValue<3>::Vector > VecField;
    typedef FieldConstant<3, FieldValue<2>::TensorFixed > TensorField;

	//Field<3, ScalarField > scalar_field;

    /*class EqData : public FieldSet {
        EqData() : FieldSet()
        {}
    };*/

	void SetUp() {
	    Profiler::initialize();

	    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

        FilePath mesh_file("mesh/simplest_cube.msh", FilePath::input_file);
        mesh_ = new Mesh;
        ifstream in(string( mesh_file ).c_str());
        mesh_->read_gmsh_from_stream(in);
	}

	void TearDown() {
		delete mesh_;
	}

	template<class Type>
	void read_input_field(string type_name) {
        Input::Type::Record region_rec("Region", "Definition of region of elements.");
        region_rec.declare_key("rid", Input::Type::Integer(), Input::Type::Default::obligatory(),"" );
        region_rec.declare_key("init_pressure", Type::input_type, Input::Type::Default::obligatory(),"" );
        region_rec.finish();

        Input::Type::Record rec_type("Test","");
        rec_type.declare_key(type_name, IT::Array(region_rec), Input::Type::Default::obligatory(),"" );
        rec_type.finish();

        Input::JSONToStorage reader( input, rec_type );
        rec_ = reader.get_root_interface<Input::Record>();
	}

	Input::Record rec_;
	Mesh *mesh_;
    Space<3>::Point point;

    FieldSet data_;

	/*double (Point&, ElAccessor&) *fce_;

	double value(Point &p, ElAccessor &elm) {
        return fce_[elm.region_idx()].fce_(p,elm);
    }*/

};


TEST_F(FieldSpeed, scalar) {
	this->read_input_field<ScalarField>("scalar_data");

	Input::Array regions = rec_.val<Input::Array>("scalar_data");

    data_.set_mesh(*mesh_);
    data_.set_input_list( regions );

    /*data_+=scalar_field
                .name("scalar_field")
                .description("Scalar field.")
                .input_default("0.0");*/

    FOR_ELEMENTS(mesh_, elm) {
    	ElementAccessor<3> ele = (*elm).element_accessor();
        //data_[ele.region_idx().idx()]->value( point, ElementAccessor<3>(mesh_, ele.index(), ele.is_boundary() ));
    }

    /*ScalarField field_ele_pressure;
	vector<double> ele_pressure;
	ele_pressure.resize(mesh_->n_elements());
	auto ele_pressure_ptr=make_shared< ScalarField >(ele_pressure, 1);
	field_ele_pressure.set_field(mesh_->region_db().get_region_set("ALL"), ele_pressure_ptr);*/

}

