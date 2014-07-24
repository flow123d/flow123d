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


static const int loop_call_count = 10000;

template <typename T>
class FieldSpeed : public testing::Test {
public:
	typedef typename Space<3>::Point Point;
	typedef typename T::return_type ReturnType;
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

        point_ = Point("1 2 3");
        point_list_.reserve(10);
        for (int i=0; i<10; i++) point_list_.push_back( Point("1 2 3") );

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

	double call_test() {
		double sum = 0.0;
		std::vector<ReturnType> value_list(10);

		START_TIMER("single_value");
		for (int i=0; i<10*loop_call_count; i++)
			FOR_ELEMENTS(this->mesh_, ele) {
				ElementAccessor<3> elm = (*ele).element_accessor();
				sum += this->field_->value( this->point_, elm);
			}
		END_TIMER("single_value");

		START_TIMER("all_values");
		for (int i=0; i<loop_call_count; i++)
			for (int j=0; j<10; j++)
				FOR_ELEMENTS(this->mesh_, ele) {
					ElementAccessor<3> elm = (*ele).element_accessor();
					sum += this->field_->value( this->point_list_[j], elm);
				}
		END_TIMER("all_values");

		START_TIMER("value_list");
		for (int i=0; i<loop_call_count; i++)
			FOR_ELEMENTS(this->mesh_, ele) {
				ElementAccessor<3> elm = (*ele).element_accessor();
				this->field_->value_list( this->point_list_, elm, value_list);
				sum += value_list[0];
			}
		END_TIMER("value_list");
		return sum;
	}


    FceType *fce_;
    ReturnType *data_;
    ReturnType data1_;
    ReturnType data2_;
    FieldAlgorithmBase<3, FieldValue<3>::Scalar > *field_;
	Mesh *mesh_;
	Point point_;
	std::vector< Point > point_list_;

    inline ReturnType value(Point &p, ElementAccessor<3> &elm) {
    	return (this->*fce_[elm.region_idx().idx()])(p,elm);
    }

};


typedef ::testing::Types< FieldValue<3>::Scalar > TestedTypes;
TYPED_TEST_CASE(FieldSpeed, TestedTypes);


TYPED_TEST(FieldSpeed, scalar) {
	double sum_func=0.0, sum_array=0.0;


	START_TIMER("Speed test");

	START_TIMER("function");
	for (int i=0; i<10*loop_call_count; i++)
		FOR_ELEMENTS(this->mesh_, ele) {
			ElementAccessor<3> elm = (*ele).element_accessor();
			sum_func += this->value( this->point_, elm);
		}
	END_TIMER("function");

	START_TIMER("array");
	for (int i=0; i<10*loop_call_count; i++)
		FOR_ELEMENTS(this->mesh_, ele) {
			ElementAccessor<3> elm = (*ele).element_accessor();
			sum_array += this->data_[elm.region_idx().idx()];
		}
	END_TIMER("array");

	END_TIMER("Speed test");

	EXPECT_DOUBLE_EQ( sum_func, sum_array );

	Profiler::instance()->output(MPI_COMM_WORLD, cout);

}

TYPED_TEST(FieldSpeed, field_constant) {
	this->field_ = new FieldConstant<3, FieldValue<3>::Scalar>();
	((FieldConstant<3, FieldValue<3>::Scalar> *)this->field_)->set_value(1.5);

	double sum_const = this->call_test();

	EXPECT_DOUBLE_EQ( sum_const, (21 * 1.5 * this->mesh_->n_elements() * loop_call_count) );

	Profiler::instance()->output(MPI_COMM_WORLD, cout);

}

