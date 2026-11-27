#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <vector>

#include "fem/eigen_tools.hh"
#include "fields/outer_prod_vec.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"


template<unsigned int dim, unsigned int spacedim>
class ShapeGradientManualHandler {
public:
    /// Constructor, set number of elements number of points on element, number of DOFs
	ShapeGradientManualHandler()
    : n_points_(4), n_dofs_(2) {
        uint n_nozero_rows = dim*n_dofs_;
        uint n_rows = n_nozero_rows + dim*n_points_ + n_points_*n_dofs_;

        data_table_.resize(n_rows);
        eigen_tools::resize_table(data_table_, spacedim);

        // fill non-zero rows
        for (uint i_row=0; i_row<n_nozero_rows; ++i_row) {
    	    for (uint i_col=0; i_col<spacedim; ++i_col) {
    	    	uint r = (i_row + i_col + 1) - (i_row / 2);
    	    	data_table_(i_row)(i_col) = 1.0*r;
            }
        }
        // fill empty (zero) rows
        for (uint i_row=n_nozero_rows; i_row<n_rows; ++i_row) {
            for (uint i_col=0; i_col<spacedim; ++i_col) {
            	data_table_(i_row)(i_col) = 0.0;
            }
        }

        // fill ref_grad_scalar_shape
        ref_grad_scalar_shape_.resize(spacedim);
        for (uint i=0; i<spacedim; ++i)
        	ref_grad_scalar_shape_[i].resize(n_points_);
        for (uint i_grad=0; i_grad<spacedim; ++i_grad) {
            for (uint i_c=0; i_c<n_points_; ++i_c) {
                ref_grad_scalar_shape_[i_grad](i_c) = 1 + 2*i_grad + i_c;
            }
        }
    }

    /// Evaluate gradients
    double eval() {
        uint grad_expand_begin = dim*n_dofs_;
        uint result_begin = grad_expand_begin + dim*n_points_;
        uint n_rows = result_begin + n_dofs_*n_points_;

        // copy ref gradient data
        for (uint i_grad=0; i_grad<dim; ++i_grad)
            for (uint i_c=0; i_c<n_points_; ++i_c) {
                ArrayDbl &data_table_row = data_table_(grad_expand_begin + i_grad*n_points_ + i_c);
                for (uint i_col=0; i_col<spacedim; ++i_col)
                    data_table_row(i_col) = ref_grad_scalar_shape_[i_grad][i_c];

            }

        // reset result rows
        for (uint i_row=result_begin; i_row<n_rows; ++i_row) {
            for (uint i_col=0; i_col<spacedim; ++i_col) {
            	data_table_(i_row)(i_col) = 0.0;
            }
        }

        // evaluate gradients over all DOFs
        auto shape_grad_value = Eigen::Map<Eigen::Matrix<ArrayDbl, 2, 4>>(data_table_.data() + result_begin, 2, 4); // result
        for (uint i=0; i<dim; ++i) {
            auto vec1 = Eigen::Map<Eigen::Matrix<ArrayDbl, 2, 1>>(data_table_.data() + 2*i, 2, 1);
            auto vec2 = Eigen::Map<Eigen::Matrix<ArrayDbl, 4, 1>>(data_table_.data() + grad_expand_begin + 4*i, 4, 1);
            shape_grad_value += vec1 * vec2.transpose();
        }

        return shape_grad_value(0,0)(0) + shape_grad_value(0,0)(1) + shape_grad_value(0,0)(2);
    }

//    /// Perform data_table_ to standard output
//    void print_table() {
//        uint n_rows = dim*n_dofs_ + dim*n_points_ + n_points_*n_dofs_;
//
//        std::cout << "Data: " << std::endl;
//        for (uint i_row=0; i_row<spacedim; ++i_row) {
//            for (uint i_col=0; i_col<n_rows; ++i_col)
//                std::cout << data_table_(i_col)(i_row) << " ";
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//    }

//    /// Perform ref_grad_scalar_shape_ to standard output
//    void print_ref_grads() {
//        std::cout << "Ref grad: " << std::endl;
//        for (uint i=0; i<spacedim; ++i) {
//            std::cout << ref_grad_scalar_shape_[i].transpose() << std::endl;
//        }
//        std::cout << std::endl;
//    }

    TableDbl data_table_;
    std::vector<Eigen::Vector<double, Eigen::Dynamic>> ref_grad_scalar_shape_;
    uint n_points_;
    uint n_dofs_;
};


template<unsigned int dim, unsigned int spacedim>
class ShapeGradientOuterProdHandler {
public:
    /// Constructor, set number of elements number of points on element, number of DOFs
	ShapeGradientOuterProdHandler() {
		// Initialize arrays
	    for (uint i=0; i<spacedim; ++i) {
	    	for (uint j=0; j<dim; ++j) {
	            arr1_(i,j).data.resize(2, 1);
	            arr1_(i,j).data << (i+j+1), (i+j+2);
	    	}
	        arr2_(i).data.resize(4, 1);
	        arr2_(i).data << (2*i+1), (2*i+2), (2*i+3), (2*i+4);
	    }

    }

    /// Evaluate gradients
    double eval() {
        result_ = arr1_ * arr2_;
        return result_(0).data(0,0) + result_(1).data(0,0) + result_(2).data(0,0);
    }

    void profiler_output(std::string file_name) {
		FilePath fp(file_name + "_profiler.json", FilePath::output_file);
		Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
	}


    Eigen::Matrix<OuterProdVec, spacedim, dim> arr1_;
    Eigen::Matrix<OuterProdVec, dim, 1> arr2_;
    Eigen::Matrix<OuterProdVec, dim, 1> result_;
};


TEST(OuterProduct, speed_test) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Profiler::instance();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    static const uint n_loops = 1e05;
    double m_result=0.0, op_result=0.0;

    ShapeGradientManualHandler<3,3> manual;
    START_TIMER("manual");
    for (uint i=0; i<n_loops; ++i)
        m_result += manual.eval();
    END_TIMER("manual");

    ShapeGradientOuterProdHandler<3,3> outer_product;
    START_TIMER("outer_product");
    for (uint i=0; i<n_loops; ++i)
        op_result += outer_product.eval();
    END_TIMER("outer_product");

    EXPECT_DOUBLE_EQ( m_result, op_result );
    outer_product.profiler_output("outer_product");
}


TEST(OuterProduct, outer_product) {
	// https://stackoverflow.com/questions/40829887/how-do-i-do-outer-product-of-tensors-in-eigen
    test_add_operator_simple();
    std::cout << "==============" << std::endl;
    test_add_operator();
    std::cout << "==============" << std::endl;
    test_multi_operator_simple();
    std::cout << "==============" << std::endl;
    test_multi_operator_vec_scalar();
    std::cout << "==============" << std::endl;
    test_multi_operator_vec_vec();
    std::cout << "==============" << std::endl;
    test_multi_operator_mat_vec();
}


