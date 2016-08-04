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
bool expect_arma_eqal(const ArmaMat1 &ref_arma_mat, const ArmaMat2 &arma_mat)
{
    auto ref_shape = mat_shape(ref_arma_mat);
    auto shape = mat_shape(arma_mat);
    if (ref_shape[0] != shape[0] || ref_shape[1] != shape[1]) {
        EXPECT_EQ(ref_shape[0], shape[0]); // rows
        EXPECT_EQ(ref_shape[1], shape[1]); // cols
        return false;
    }
    double magnitude = std::max( arma::norm(ref_arma_mat, 1), arma::norm(arma_mat, 1) );
    double error = arma::norm(ref_arma_mat - arma_mat, 1)/magnitude;
    if (error > 8*std::numeric_limits<double>::epsilon()) {
        unsigned int w = 10* arma_mat.n_cols;
        std::cerr << std::setw(w) << "Expected" << std::setw(w) << "Result" << std::endl;
        for(unsigned int i_row = 0; i_row < arma_mat.n_rows; i_row++) {
            std::cerr << std::setw(10) << ref_arma_mat.row(i_row) << std::setw(10) << arma_mat.row(i_row) << std::endl;
        }
        return false;
    }
    return true;
}

#define EXPECT_ARMA_EQ( A, B) \
    EXPECT_TRUE(expect_arma_eqal(A, B))



#endif /* UNIT_TESTS_ARMA_EXPECT_HH_ */
