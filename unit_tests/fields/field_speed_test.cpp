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


using namespace std;


class FieldSpeed : public testing::Test {
public:
	typedef typename Space<3>::Point Point;
	typedef FieldValue<3>::Scalar::return_type ReturnType;
    typedef ReturnType(FieldSpeed::*FceType)(Point&, ElementAccessor<3>&);

    ReturnType fce1(Point &p, ElementAccessor<3> &elm) {
    	return data1_;
    }

    ReturnType fce2(Point &p, ElementAccessor<3> &elm) {
    	return data2_;
    }

	void SetUp() {
	    Profiler::initialize();

	    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

        FilePath mesh_file("mesh/simplest_cube.msh", FilePath::input_file);
        mesh_ = new Mesh;
        ifstream in(string( mesh_file ).c_str());
        mesh_->read_gmsh_from_stream(in);

        point = Point("1 2 3");

    	fce_ = new FceType[mesh_->region_db().size()];
    	data_ = new ReturnType[mesh_->region_db().size()];

        data1_ = 1.1;
    	data2_ = 2.2;

    	data_[0] = data1_;
    	data_[1] = data2_;
    	data_[2] = data2_;
    	data_[3] = data1_;
    	data_[4] = data2_;
    	data_[5] = data1_;
    	data_[6] = data1_;
    	data_[7] = data2_;

    	fce_[0] = (&FieldSpeed::fce1);
        fce_[1] = (&FieldSpeed::fce2);
        fce_[2] = (&FieldSpeed::fce2);
        fce_[3] = (&FieldSpeed::fce1);
        fce_[4] = (&FieldSpeed::fce2);
        fce_[5] = (&FieldSpeed::fce1);
        fce_[6] = (&FieldSpeed::fce1);
        fce_[7] = (&FieldSpeed::fce2);
	}

	void TearDown() {
		Profiler::uninitialize();

		delete mesh_;
	}


	static const int loop_call_count = 100000;
    FceType *fce_;
    ReturnType *data_;
    ReturnType data1_;
    ReturnType data2_;
    FieldAlgorithmBase<3, FieldValue<3>::Scalar > *field_;
	Mesh *mesh_;
	Point point;

    inline ReturnType value(Point &p, ElementAccessor<3> &elm) {
    	return (this->*fce_[elm.region_idx().idx()])(p,elm);
    }

};


TEST_F(FieldSpeed, scalar) {
	double sum_func=0.0, sum_array=0.0, sum_const=0.0;

	START_TIMER("Speed test");

	START_TIMER("function");
	for (int i=0; i<loop_call_count; i++)
		FOR_ELEMENTS(mesh_, ele) {
			ElementAccessor<3> elm = (*ele).element_accessor();
			sum_func += this->value( point, elm);
		}
	END_TIMER("function");

	START_TIMER("array");
	for (int i=0; i<loop_call_count; i++)
		FOR_ELEMENTS(mesh_, ele) {
			ElementAccessor<3> elm = (*ele).element_accessor();
			sum_array += this->data_[elm.region_idx().idx()];
		}
	END_TIMER("array");

	START_TIMER("field_const");
	field_ = new FieldConstant<3, FieldValue<3>::Scalar>();
	((FieldConstant<3, FieldValue<3>::Scalar> *)field_)->set_value(1.5);
	for (int i=0; i<loop_call_count; i++)
		FOR_ELEMENTS(mesh_, ele) {
			ElementAccessor<3> elm = (*ele).element_accessor();
			sum_const += field_->value( point, elm);
		}
	END_TIMER("field_const");

	END_TIMER("Speed test");

	EXPECT_DOUBLE_EQ( sum_func, sum_array );
	EXPECT_DOUBLE_EQ( sum_const, (1.5 * mesh_->n_elements() * loop_call_count) );

	Profiler::instance()->output(MPI_COMM_WORLD, cout);

}

