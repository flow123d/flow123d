#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include <vector>

#include "fem/eigen_tools.hh"
#include "fields/outer_prod_vec.hh"


template<unsigned int dim, unsigned int spacedim>
class PatchData {
public:
    /// Constructor, set number of elements number of points on element, number of DOFs
    PatchData(uint n_elem, uint n_elem_points, uint n_dofs)
    : n_elem_(n_elem), n_elem_points_(n_elem_points), n_dofs_(n_dofs) {
        uint n_nozero_rows = dim*spacedim;
        uint n_rows = n_nozero_rows + dim*n_dofs_ + spacedim*n_dofs_;
        uint n_points = n_elem_ * n_elem_points_;

        data_table_.resize(n_rows);
        eigen_tools::resize_table(data_table_, n_points);

        // fill rows simulate inverse jacobian
        for (uint i_row=0; i_row<n_nozero_rows; ++i_row) {
    	    for (uint i_pt=0; i_pt<n_points; ++i_pt) {
    	    	uint r = (10 + 3*i_row*i_pt + 2*i_row - i_pt) % 10;
    	    	data_table_(i_row)(i_pt) = 0.11 + 0.1*r;
            }
        }
        // fill empty rows (result)
        for (uint i_row=n_nozero_rows; i_row<n_rows; ++i_row) {
            for (uint i_pt=0; i_pt<n_points; ++i_pt) {
            	data_table_(i_row)(i_pt) = 0.0;
            }
        }

        // fill ref_grad_scalar_shape
        ref_grad_scalar_shape_.resize(n_elem_points_);
        for (uint i=0; i<n_elem_points_; ++i)
        	ref_grad_scalar_shape_[i].resize(n_dofs_);
        for (uint i_pt=0; i_pt<n_elem_points_; ++i_pt) {
            for (uint i_dof=0; i_dof<n_dofs_; ++i_dof) {
                uint r = (10 + 3*i_pt*i_dof + 2*i_pt - i_dof) % 10;
                for (uint i_c=0; i_c<dim; ++i_c)
                    ref_grad_scalar_shape_[i_pt][i_dof](i_c) = 0.13 + 0.02*i_c + 0.1*r;
            }
        }
    }

    /// Perform data_table_ to standard output
    void print_table() {
        uint n_rows = dim*spacedim + dim*n_dofs_ + spacedim*n_dofs_;
        uint n_points = n_elem_ * n_elem_points_;

        std::cout << "Data: " << std::endl;
        for (uint i_row=0; i_row<n_points; ++i_row) {
            for (uint i_col=0; i_col<n_rows; ++i_col)
                std::cout << data_table_(i_col)(i_row) << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    /// Perform ref_grad_scalar_shape_ to standard output
    void print_ref_grads() {
        std::cout << "Ref grad: " << std::endl;
        for (uint i_pt=0; i_pt<n_elem_points_; ++i_pt) {
            for (uint i_dof=0; i_dof<n_dofs_; ++i_dof)
                std::cout << "(" << ref_grad_scalar_shape_[i_pt][i_dof](0) << ", " << ref_grad_scalar_shape_[i_pt][i_dof](1) << ", " << ref_grad_scalar_shape_[i_pt][i_dof](2) << ") ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    /// Evaluate gradients
    void eval() {
        uint inv_jac_begin = 0;
        uint grad_expand_begin = inv_jac_begin + dim*spacedim;
        uint result_begin = grad_expand_begin + dim*n_dofs_;

        // copy ref gradient data
        for (uint i_dof=0; i_dof<n_dofs_; ++i_dof)
            for (uint i_c=0; i_c<dim; ++i_c) {
                ArrayDbl &data_table_row = data_table_(grad_expand_begin + i_dof*dim+  i_c);
                for (uint i_pt=0; i_pt<data_table_row.rows(); ++i_pt)
                    data_table_row(i_pt) = ref_grad_scalar_shape_[i_pt % n_elem_points_][i_dof][i_c];
            }

        // evaluate gradients over all DOFs
        auto inv_jac_value = Eigen::Map<Eigen::Matrix<ArrayDbl, dim, spacedim>>(data_table_.data() + inv_jac_begin, dim, spacedim);
        for (uint i_dof=0; i_dof<n_dofs_; ++i_dof) {
            auto ref_shape_grad_value = Eigen::Map<Eigen::Matrix<ArrayDbl, dim, 1>>(data_table_.data() + grad_expand_begin + dim*i_dof, dim, 1);
            auto shape_grad_value = Eigen::Map<Eigen::Matrix<ArrayDbl, spacedim, 1>>(data_table_.data() + result_begin + spacedim*i_dof, spacedim, 1); // result
            shape_grad_value = inv_jac_value.transpose() * ref_shape_grad_value;
        }
    }

    TableDbl data_table_;
    std::vector< std::vector<arma::vec3> > ref_grad_scalar_shape_;
    uint n_elem_;
    uint n_elem_points_;
    uint n_dofs_;
};

TEST(OuterProduct, manual_result) {
    PatchData<3,3> pd(2, 4, 2);
    pd.print_table();
    pd.print_ref_grads();
    pd.eval();
    pd.print_table();
}


//TEST(OuterProduct, outer_product_fast) {
//    Eigen::Array<Eigen::Index, 0, 0> empty_index_list = {};
//    Eigen::Matrix<double, 3, 3> A_ij;
//    Eigen::Matrix<double, 3, 1> B_kl;
//    auto C_ijkl = A_ij.contract(B_kl, empty_index_list);
//}


TEST(OuterProduct, outer_product) {
	// https://stackoverflow.com/questions/40829887/how-do-i-do-outer-product-of-tensors-in-eigen
    test_add_operator_simple();
    std::cout << "==============" << std::endl;
    test_multi_operator_simple();
    std::cout << "==============" << std::endl;
    test_add_operator();
    std::cout << "==============" << std::endl;
    test_multi_operator();
    std::cout << "==============" << std::endl;
    test_multi_mat_vec_operator();
}


