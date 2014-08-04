/*
 * field_speed_test.cpp
 *
 *  Created on: Feb 3, 2013
 *      Author: jb
 */

#define TEST_USE_MPI

#include <flow_gtest_mpi.hh>

#include "fields/field_constant.hh"
#include "fields/field_formula.hh"
#include "fields/field_python.hh"
#include "fields/field_elementwise.hh"
#include "fields/field.hh"
#include "fields/field_set.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"


using namespace std;


static const int loop_call_count = 1;
static const int list_size = 10;


template <typename T>
class FieldSpeed : public testing::Test {
public:
	typedef typename Space<3>::Point Point;
	typedef FieldAlgorithmBase<3, T > FieldBaseType;
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

        set_values();
	}

	void TearDown() {
		Profiler::uninitialize();

		delete mesh_;
	}

	ReturnType call_test() {


		START_TIMER("single_value");
		for (int i=0; i<list_size*loop_call_count; i++)
			FOR_ELEMENTS(this->mesh_, ele) {
				ElementAccessor<3> elm = (*ele).element_accessor();
				test_result_sum_ += this->field_->value( this->point_, elm);
			}
		END_TIMER("single_value");

		START_TIMER("all_values");
		for (int i=0; i<loop_call_count; i++)
			for (int j=0; j<list_size; j++)
				FOR_ELEMENTS(this->mesh_, ele) {
					ElementAccessor<3> elm = (*ele).element_accessor();
					test_result_sum_ += this->field_->value( this->point_list_[j], elm);
				}
		END_TIMER("all_values");

		START_TIMER("value_list");
		for (int i=0; i<loop_call_count; i++)
			FOR_ELEMENTS(this->mesh_, ele) {
				ElementAccessor<3> elm = (*ele).element_accessor();
				this->field_->value_list( this->point_list_, elm, value_list);
				test_result_sum_ += value_list[0];
			}
		END_TIMER("value_list");
		return test_result_sum_;
	}

	void set_values() {
        n_comp_ = 3;

        point_ = Point("1 2 3");
        point_list_.reserve(10);
        for (int i=0; i<10; i++) point_list_.push_back( Point("1 2 3") );

    	fce_ = new FceType[mesh_->region_db().size()];
    	data_ = new ReturnType[mesh_->region_db().size()];

    	set_data(data1_);

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

	void set_data(FieldValue<3>::Scalar::return_type val) {
        data1_ = 1.25;
    	data2_ = 2.50;
    	field_const_val_ = 1.25;
    	field_elementwise_val_ = 4.5;
    	test_result_sum_ = 0.0;
    	input_type_name_ = "scalar";
    	value_list= std::vector<ReturnType>(list_size);
	}

	void set_data(FieldValue<3>::Vector::return_type val) {
		data1_ = arma::vec("1.25 2.25 3.25");
		data2_ = arma::vec("2.50 4.50 6.50");
		field_const_val_ = arma::vec("1.75 2.75 3.75");
		field_elementwise_val_ = arma::vec("9 18 27");
		test_result_sum_ = arma::vec("0.0 0.0 0.0");
		input_type_name_ = "vector";
		value_list= std::vector<ReturnType>(list_size, ReturnType(n_comp_,1));
	}

	void set_data(FieldValue<3>::VectorFixed::return_type val) {
		data1_ = arma::vec3("1.25 3.75 6.25");
		data2_ = arma::vec3("2.50 7.50 12.50");
		field_const_val_ = arma::vec3("1.75 3.75 5.75");
		field_elementwise_val_ = arma::vec3("9 18 27");
		test_result_sum_ = arma::vec3("0.0 0.0 0.0");
		input_type_name_ = "vector_fixed";
		value_list= std::vector<ReturnType>(list_size);
	}

	void test_result(FieldValue<3>::Scalar::return_type expected, double multiplicator) {
		cout << "Test result: " << this->test_result_sum_ << endl;
		cout << "Expected:    " << multiplicator * expected << endl;
		EXPECT_DOUBLE_EQ( this->test_result_sum_, multiplicator * expected * loop_call_count );
	}

	void test_result(FieldValue<3>::Vector::return_type expected, double multiplicator) {
		cout << "Test result: " << this->test_result_sum_ << endl;
		cout << "Expected:    " << multiplicator * expected << endl;
		for (int i=0; i<3; i++) {
			EXPECT_DOUBLE_EQ( this->test_result_sum_[i], multiplicator * expected[i] * loop_call_count );
		}
	}

	void test_result(FieldValue<3>::VectorFixed::return_type expected, double multiplicator) {
		cout << "Test result: " << this->test_result_sum_ << endl;
		cout << "Expected:    " << multiplicator * expected << endl;
		for (int i=0; i<3; i++) {
			EXPECT_DOUBLE_EQ( this->test_result_sum_[i], multiplicator * expected[i] * loop_call_count );
		}
	}

	void profiler_output() {
		static ofstream os( FilePath("speed_test_profiler.log", FilePath::output_file) );
		Profiler::instance()->output(MPI_COMM_WORLD, os);
		os << "" << std::setfill('=') << setw(80) << "" << std::setfill(' ') << endl << endl;
	}

	void read_input(const string &input, const string &field_name) {
	    /*Input::Type::Record  rec_type("FieldSpeedTest","");
	    rec_type.declare_key(field_name, Input::Type::Array(FieldBaseType::input_type), Input::Type::Default::obligatory(),"" );
	    rec_type.finish();

	    //Input::Type::Array array_type = Input::Type::Array(FieldBaseType::input_type);

	    Input::JSONToStorage reader( input, rec_type );
	    Input::Record in_rec=reader.get_root_interface<Input::Record>();

	    this->field_.name("data");
	    this->field_.set_mesh(*(this->mesh_));
	    this->field_.set_input_list( in_rec.val<Input::Array>(field_name) );
	    TimeGovernor tg(0.0, 1.0);
	    this->field_.set_time(tg);*/
	}


    FceType *fce_;
    ReturnType *data_;
    ReturnType data1_;
    ReturnType data2_;
    ReturnType test_result_sum_;
    ReturnType field_const_val_;
    ReturnType field_elementwise_val_;
    std::vector<ReturnType> value_list;
    FieldAlgorithmBase<3, T> *field_;
    //Field<3, T> field_;
	Mesh *mesh_;
	Point point_;
	std::vector< Point > point_list_;
	string input_type_name_;
	unsigned int n_comp_;

    inline ReturnType value(Point &p, ElementAccessor<3> &elm) {
    	return (this->*fce_[elm.region_idx().idx()])(p,elm);
    }

};


typedef ::testing::Types< FieldValue<3>::Scalar, FieldValue<3>::Vector, FieldValue<3>::VectorFixed > TestedTypes;
//typedef ::testing::Types< FieldValue<3>::Vector > TestedTypes;
TYPED_TEST_CASE(FieldSpeed, TestedTypes);

TYPED_TEST(FieldSpeed, array) {
    START_TIMER("array");
	START_TIMER("single_value");
	for (int i=0; i<10*loop_call_count; i++)
		FOR_ELEMENTS(this->mesh_, ele) {
			ElementAccessor<3> elm = (*ele).element_accessor();
			this->test_result_sum_ += this->data_[elm.region_idx().idx()];
		}
	END_TIMER("single_value");
	END_TIMER("array");

	this->test_result(this->data1_, 130);
	this->profiler_output();
}


TYPED_TEST(FieldSpeed, virtual_function) {
	START_TIMER("virtual_function");
	START_TIMER("single_value");
	for (int i=0; i<10*loop_call_count; i++)
		FOR_ELEMENTS(this->mesh_, ele) {
			ElementAccessor<3> elm = (*ele).element_accessor();
			this->test_result_sum_ += this->value( this->point_, elm);
		}
	END_TIMER("single_value");

	START_TIMER("all_values");
	for (int i=0; i<loop_call_count; i++)
		for (int j=0; j<10; j++)
			FOR_ELEMENTS(this->mesh_, ele) {
				ElementAccessor<3> elm = (*ele).element_accessor();
				this->test_result_sum_ += this->value( this->point_list_[j], elm);
			}
	END_TIMER("all_values");
	END_TIMER("virtual_function");

	this->test_result(this->data1_, 2 * 130);
	this->profiler_output();
}



string constant_input = R"JSON(
{   
    constant_expr_scalar=[
        { TYPE="FieldConstant", value=1.25 }
    ],
    constant_expr_vector=[
        { TYPE="FieldConstant", value=[1.75, 2.75, 3.75] }
    ],
    constant_expr_vector_fixed=[
        { TYPE="FieldConstant", value=[1.75, 3.75, 5.75] }
    ]
}
)JSON";

TYPED_TEST(FieldSpeed, field_constant) {
	//string key_name = "constant_expr_" + this->input_type_name_;
	//this->read_input(constant_input, key_name);

    this->field_ = new FieldConstant<3, TypeParam>(this->n_comp_);
	((FieldConstant<3, TypeParam> *)this->field_)->set_value(this->field_const_val_);

	START_TIMER("field_constant");
	this->call_test();
	END_TIMER("field_constant");

	this->test_result( this->field_const_val_, (21 * this->mesh_->n_elements()) );
	this->profiler_output();
}


string formula_input = R"INPUT(
{   
   constant_expr_scalar={ 
       TYPE="FieldFormula",
       value="1.25"
   },
   constant_expr_vector={ 
       TYPE="FieldFormula",
       value=["1.75", "2.75", "3.75"]
   },
   constant_expr_vector_fixed={ 
       TYPE="FieldFormula",
       value=["1.75", "3.75", "5.75"]
   },
   simple_expr={ 
       TYPE="FieldFormula",
       value="x^2"
   },
   full_expr={ 
       TYPE="FieldFormula",
       value="x+y+z+x^2+y^2+z^2"
   }
}
)INPUT";

TYPED_TEST(FieldSpeed, field_formula) {
	typedef FieldAlgorithmBase<3, TypeParam > FieldBaseType;
	string key_name = "constant_expr_" + this->input_type_name_;

    Input::Type::Record rec_type("FieldFormulaTest","");
    rec_type.declare_key(key_name, FieldBaseType::input_type, Input::Type::Default::obligatory(), "" );
    rec_type.finish();

    // read input string
    Input::JSONToStorage reader( formula_input, rec_type );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

	this->field_ = new FieldFormula<3, TypeParam>(this->n_comp_);
	this->field_->init_from_input(in_rec.val<Input::AbstractRecord>(key_name));
	this->field_->set_time(0.0);

	START_TIMER("field_formula");
	this->call_test();
	END_TIMER("field_formula");

	cout << "Test result: " << this->test_result_sum_ << endl;
	//this->test_result( this->field_const_val_, (9 * 21) );
	this->profiler_output();
}


#ifdef HAVE_PYTHON
string python_input = R"CODE(
def func_scalar(x,y,z):
    return ( 1.25, )     # one value tuple

def func_vector(x,y,z):
    return ( 1.75, 2.75, 3.75 )     # one value tuple

def func_vector_fixed(x,y,z):
    return ( 1.75, 3.75, 5.75 )     # one value tuple

)CODE";

TYPED_TEST(FieldSpeed, field_python) {
    this->field_ = new FieldPython<3, TypeParam>(this->n_comp_);
    ((FieldPython<3, TypeParam> *)this->field_)->set_python_field_from_string(python_input, "func_" + this->input_type_name_);
    this->field_->set_time(0.0);

	START_TIMER("field_python");
	this->call_test();
	END_TIMER("field_python");

	this->test_result( this->field_const_val_, (9 * 21) );
	this->profiler_output();
}
#endif // HAVE_PYTHON

string elementwise_input = R"INPUT(
{   
   scalar={
       TYPE="FieldElementwise",
       gmsh_file="fields/simplest_cube_data.msh",
       field_name="scalar"
   },
   vector={
       TYPE="FieldElementwise",
       gmsh_file="fields/simplest_cube_data.msh",
       field_name="vector_fixed"
   },
   vector_fixed={
       TYPE="FieldElementwise",
       gmsh_file="fields/simplest_cube_data.msh",
       field_name="vector_fixed"
   }
}
)INPUT";

TYPED_TEST(FieldSpeed, field_elementwise) {
	Input::Type::Record  rec_type("Test","");
    rec_type.declare_key("scalar", FieldElementwise<3, FieldValue<3>::Scalar >::input_type, Input::Type::Default::obligatory(),"" );
    rec_type.declare_key("vector_fixed", FieldElementwise<3, FieldValue<3>::VectorFixed >::input_type, Input::Type::Default::obligatory(),"" );
    rec_type.declare_key("vector", FieldElementwise<3, FieldValue<3>::Vector >::input_type, Input::Type::Default::obligatory(),"" );
    rec_type.finish();

    Input::JSONToStorage reader( elementwise_input, rec_type );
    Input::Record rec=reader.get_root_interface<Input::Record>();

    this->field_ = new FieldElementwise<3, TypeParam>(this->n_comp_);
    this->field_->init_from_input( rec.val<Input::Record>(this->input_type_name_) );
    this->field_->set_mesh(this->mesh_,false);
    this->field_->set_time(0.0);

    START_TIMER("field_elementwise");
	this->call_test();
	END_TIMER("field_elementwise");

	this->test_result( this->field_elementwise_val_, 21 );
	this->profiler_output();
}
