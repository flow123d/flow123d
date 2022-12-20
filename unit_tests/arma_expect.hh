/*
 * arma_expect.hh
 *
 *  Created on: Jul 4, 2016
 *      Author: jb
 */

#ifndef UNIT_TESTS_ARMA_EXPECT_HH_
#define UNIT_TESTS_ARMA_EXPECT_HH_

#include <armadillo>
#include <numeric>
#include <iostream>
#include <algorithm>

template <int N, int M>
arma::uvec mat_shape(const arma::mat::fixed<N, M> &x) {
    return arma::uvec({ N, M});
}

template <class Mat>
arma::uvec mat_shape(const Mat &x) {
    return arma::uvec({ x.n_rows, x.n_cols});
}



template<class ArmaMat1, class ArmaMat2>
bool _expect_arma_eqal(const ArmaMat1 &ref_arma_mat, const ArmaMat2 &arma_mat, std::ostream &fail_message)
{
    bool no_failure=true;
    auto ref_shape = mat_shape(ref_arma_mat);
    auto shape = mat_shape(arma_mat);
    if (ref_shape[0] != shape[0] || ref_shape[1] != shape[1]) {
        if (ref_shape[0] != shape[0]) fail_message << "Different number of rows of field value and ref value!" << std::endl;
        if (ref_shape[1] != shape[1]) fail_message << "Different number of cols of field value and ref value!" << std::endl;
        no_failure=false;
    } else {
        double magnitude = std::max( arma::norm(ref_arma_mat, 1), arma::norm(arma_mat, 1) );
        // abs criterium
        if (magnitude < 8*std::numeric_limits<double>::epsilon()) return true;

        double error = arma::norm(ref_arma_mat - arma_mat, 1)/magnitude;
        // rel criterium
        if (error > 8*std::numeric_limits<double>::epsilon()) {
            unsigned int w = 11* arma_mat.n_cols;
            fail_message << std::setw(w) << "Expected" << std::setw(w) << "Result"
                    << " rel. error: " << error << std::endl;
            for(unsigned int i_row = 0; i_row < arma_mat.n_rows; i_row++) {
                fail_message << std::setw(11) << ref_arma_mat.row(i_row)
                             << std::setw(11) << arma_mat.row(i_row)
                             << std::setw(20) << ref_arma_mat.row(i_row) - arma_mat.row(i_row) << std::endl;
            }
            no_failure=false;
        }
    }

    return no_failure;
}


template<class ArmaMat1, class ArmaMat2>
bool expect_arma_eqal(const ArmaMat1 &ref_arma_mat, const ArmaMat2 &arma_mat)
{
    std::ostringstream fail_message;
    bool no_failure = _expect_arma_eqal(ref_arma_mat, arma_mat, fail_message);

    if (! no_failure) {
        GTEST_NONFATAL_FAILURE_( fail_message.str().c_str() );
    }
    return no_failure;
}

#define EXPECT_ARMA_EQ( A, B ) \
    EXPECT_TRUE(expect_arma_eqal(A, B));

#define ASSERT_ARMA_EQ( A, B) \
    ASSERT_TRUE(expect_arma_eqal(A, B));


#endif /* UNIT_TESTS_ARMA_EXPECT_HH_ */
