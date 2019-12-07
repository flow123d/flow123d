#include <flow_gtest.hh>
#include "system/armor.hh"
#include "system/logger.hh"
#include <armadillo>

/**
 * TODO:
 * - Compiler can  not find operators since 'MatInstance...' is not derived from 'const Mat...' ??
 * - armadillo defines fixed inside of Mat, is this a reason? Try it:
 *   class Mat<T, nr, nc>::Data;
 * - is MatInstance what we need? Since we always use arma::Mat as temporary value (results of operations)
 * - make separate file for experiments
 *
 */
//template<class Type, int nRows, int nCols>
//class MatInstance : public Armor::Mat<Type, nRows, nCols> {
//    typedef typename Armor::Mat<Type, nRows, nCols>::Arma Arma;
//    Type data[16];
//public:
//    inline MatInstance()
//    : Armor::Mat<Type, nRows, nCols>(data)
//    {}
//
//    inline MatInstance(std::initializer_list<std::initializer_list<Type>> list)
//    // construct from list of columns
//    : MatInstance()
//    {
//        Armor::Mat<Type, nRows, nCols>::operator=(Arma(list));
//    }
//
//    inline MatInstance(std::initializer_list<Type> list)
//    // construct from list of all entries, the fastest increase in rows
//    : MatInstance()
//    {
//        Armor::Mat<Type, nRows, nCols>::operator=(Arma(list));
//    }
//
//
//    inline const MatInstance<Type, nRows, nCols> & operator=(const MatInstance<Type, nRows, nCols> & other)
//    {
//        Armor::Mat<Type, nRows, nCols>::operator=(other);
//        return *this;
//    }
//
//    inline const MatInstance<Type, nRows, nCols> & operator=(const Arma & other) {
//        Armor::Mat<Type, nRows, nCols>::operator=(other);
//        return *this;
//    }
//};

template<class Type, int nr, int nc>
Armor::Mat<Type, nr, nc> make_mat(std::initializer_list<std::initializer_list<Type>> list) {
    //DebugOut() << "size: " << nr * nc << "\n";
    //auto ptr = std::make_shared<Type>(2 * nr * nc);
    Type * ptr = new Type[nr * nc];
    //Armor::Mat<Type, nr, nc>  mat(ptr.get());
    Armor::Mat<Type, nr, nc>  mat(ptr);
    //DebugOut() << "Size: " << mat.n_rows * mat.n_cols << "\n";
    mat = typename Armor::Mat<Type, nr, nc>::Arma(list);
    return mat;
}


/**************** Test MatInstance **************/

TEST(Armor_test, direct_construction) {
    typename Armor::Mat<double, 3, 1>::Arma a31({1, 2, 3});
    Armor::Mat<double, 3, 1> m31(a31);
    Armor::Mat<double, 3, 1> m32(a31);
    m31 = m32;
    m31 = a31;
    EXPECT_EQ(1, m31[0]);
    EXPECT_EQ(2, m31[1]);
    EXPECT_EQ(3, m31[2]);
}

TEST(Armor_test, constructor_list1) {
    auto m1 = make_mat<double, 3, 1>({{1}, {2}, {3}});
    EXPECT_EQ(1, m1[0]);
    EXPECT_EQ(2, m1[1]);
    EXPECT_EQ(3, m1[2]);
}


TEST(Armor_test, constructor_list2) {
    arma::Mat<double>::fixed<3, 3> a1{1, 2, 3, 4, 5, 6, 7, 8, 9};
    EXPECT_EQ(4, a1.at(0, 1) );     // ROW, COLUMN
    EXPECT_EQ(2, a1.at(1, 0) );
    EXPECT_EQ(6, a1.at(2, 1) );

    arma::Mat<double>::fixed<3, 3> a2{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    EXPECT_EQ(2, a2.at(0, 1) );
    EXPECT_EQ(4, a2.at(1, 0) );
    EXPECT_EQ(8, a2.at(2, 1) );

//    auto m1 = make_mat<double, 3, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9});
//    EXPECT_EQ(4, m1.at(0, 1) );
//    EXPECT_EQ(2, m1.at(1, 0) );
//    EXPECT_EQ(6, m1.at(2, 1) );

    auto m2 = make_mat<double, 3, 3>({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    EXPECT_EQ(2, a2.at(0, 1) );
    EXPECT_EQ(4, a2.at(1, 0) );
    EXPECT_EQ(8, a2.at(2, 1) );

}


/** construction and assignment from Arma types ********/
TEST(Armor_test, from_arma_33) {
    Armor::Mat<double, 3, 3>::Arma a33 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    auto m33 = make_mat<double, 3, 3>({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    EXPECT_TRUE(Armor::is_close(
            Armor::Mat<double, 3, 3>(a33), m33));
    auto m33_none = make_mat<double, 3, 3>({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    m33_none = a33;
    EXPECT_TRUE(Armor::is_close(
            m33_none, m33));
}

TEST(Armor_test, from_arma_31) {
    Armor::Mat<double, 3, 1>::Arma a31 = {1, 2, 3};
    auto m31 = make_mat<double, 3, 1>({{1}, {2}, {3}});
    EXPECT_TRUE(Armor::is_close(
            Armor::Mat<double, 3, 1>(a31), m31));
    auto m_none = make_mat<double, 3, 1>({{1}, {2}, {3}});
    m_none = a31;
    EXPECT_TRUE(Armor::is_close(m_none, m31));

    Armor::Mat<double, 3, 1>::ArmaVec ac3 = {1, 2, 3};
    EXPECT_TRUE(Armor::is_close(
            Armor::Mat<double, 3, 1>(ac3), m31));
    m_none = make_mat<double, 3, 1>({{1}, {2}, {3}});
    m_none = ac3;
    EXPECT_TRUE(Armor::is_close(m_none, m31));
}


TEST(Armor_test, from_arma_11) {
    Armor::Mat<double, 1, 1>::Arma a11 = {3};
    auto m11 = make_mat<double, 1, 1>({{3}});
    EXPECT_TRUE(Armor::is_close(
            Armor::Mat<double, 1, 1>(a11), m11));
    auto m_none = make_mat<double, 1, 1>({{0}});
    m_none = a11;
    EXPECT_TRUE(Armor::is_close(m_none, m11));

    Armor::Mat<double, 1, 1>::ArmaVec ac1 = {3};
    EXPECT_TRUE(Armor::is_close(
            Armor::Mat<double, 1, 1>(ac1), m11));
    m_none = make_mat<double, 1, 1>({{0}});
    m_none = ac1;
    EXPECT_TRUE(Armor::is_close(m_none, m11));

    Armor::Mat<double, 1, 1>::Scalar as = 3;
    EXPECT_TRUE(Armor::is_close(
            Armor::Mat<double, 1, 1>(as), m11));
    m_none = make_mat<double, 1, 1>({{0}});
    m_none = as;
    EXPECT_TRUE(Armor::is_close(m_none, m11));

}


TEST(Armor_test, indexing) {
  // Read access.
  arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
  Armor::Mat<double, 3, 3> m1(a1);
  EXPECT_EQ(2, m1(0, 1));
  EXPECT_EQ(6, m1(1, 2));
  EXPECT_EQ(8, m1(2, 1));
  EXPECT_EQ(2, m1.at(0, 1));
  EXPECT_EQ(6, m1.at(1, 2));
  EXPECT_EQ(8, m1.at(2, 1));
  EXPECT_EQ(2, m1[3]);
  EXPECT_EQ(6, m1[7]);
  EXPECT_EQ(8, m1[5]);

  // Write access.
  m1(0, 1) = 1;
  EXPECT_EQ(1, m1(0, 1));
  EXPECT_EQ(1, m1[3]);

  m1.at(0, 1) = 3.4;
  EXPECT_EQ(3.4, m1(0, 1));
  EXPECT_EQ(3.4, m1[3]);

  m1[3] = 5.4;
  EXPECT_EQ(5.4, m1(0, 1));
  EXPECT_EQ(5.4, m1[3]);
}


/******************/

TEST(Armor_test, size) {
    auto m1 = make_mat<double, 3, 1>({{1}, {2}, {3}});
    auto m2 = make_mat<double, 3, 3>({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    EXPECT_EQ(3, m1.size());
    EXPECT_EQ(9, m2.size());
}

TEST(Armor_test, memptr) {
  arma::Mat<double>::fixed<3, 3> a1({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
  Armor::Mat<double, 3, 3> m1(a1);
  EXPECT_EQ(m1(0,0), a1(0,0));
  EXPECT_EQ(m1.memptr(), a1.memptr());
  EXPECT_EQ(&(m1(0,0)), &(a1(0,0)));
}







//TEST(Armor_test, equality_operator) {
//	MatInstance<double, 3, 1> m1;
//	arma::Mat<double>::fixed<3, 1> a1;
//	a1[0] = m1[0] = 1;
//	a1[1] = m1[1] = 2;
//	a1[2] = m1[2] = 3;
//	EXPECT_TRUE(m1 == a1);
//	a1[1] = 7;
//	EXPECT_FALSE(m1 == a1);
//}
//
//
//

//TEST(Armor_test, assignment_list1) {
//	MatInstance<double, 3, 1> m1;
//	m1 = {1, 2, 3};
//	EXPECT_EQ(1, m1[0]);
//	EXPECT_EQ(2, m1[1]);
//	EXPECT_EQ(3, m1[2]);
//}
//
//TEST(Armor_test, assignment_list2) {
//	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
//	MatInstance<double, 3, 3> m1;
//	m1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	EXPECT_TRUE(m1 == a1);
//}

//TEST(Armor_test, assignment_armor) {
//	MatInstance<double, 3, 3> m1;
//	MatInstance<double, 3, 3> m2{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	m1 = m2;
//	EXPECT_TRUE(m1 == m2);
//}
//
//TEST(Armor_test, assignment_arma) {
//	MatInstance<double, 3, 3> m1;
//	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
//	m1 = a1;
//	EXPECT_TRUE(m1 == a1);
//}
//
//
//TEST(Armor_test, arma) {
//	MatInstance<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
//    arma::Mat<double>::fixed<3, 3> a2 = m1.arma();
//    bool res = true;
//    for (uint i{0}; i < 9; ++i) {
//        if (a1[i] != a2[i]) {
//            res = false;
//        }
//    }
//	EXPECT_TRUE(res);
//}
//
//TEST(Armor_test, dot) {
//	MatInstance<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
//	double a = dot(a1, a1);
//	double m = dot(m1, m1);
//	EXPECT_EQ(a, m);
//}
//
//TEST(Armor_test, plus) {
//	MatInstance<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	Armor::Mat<double, 3, 3>::Arma a1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	Armor::Mat<double, 3, 3>::Arma res = (m1 + a1) + (a1 + m1) + (m1 + m1);
//	Armor::Mat<double, 3, 3>::Arma a6 = 6 * a1;
//	EXPECT_TRUE(arma::approx_equal(res, a6, "absdiff", 1e-6));
//}

//TEST(Armor_test, minus) {
//	MatInstance<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
//	a1 = a1 - a1;
//	m1 = m1 - m1;
//	EXPECT_TRUE(m1 == a1);
//}
//
//TEST(Armor_test, multiplication) {
//	MatInstance<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
//	a1 = a1 * a1;
//	m1 = m1 * m1;
//	EXPECT_TRUE(m1 == a1);
//}

//TEST(Armor_test, multiplication_per_elements) {
//	MatInstance<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
//	a1 = a1 % a1;
//	m1 = m1 % m1;
//	EXPECT_TRUE(m1 == a1);
//}

/*
TEST(Armor_test, array) {
	MatInstance<double, 3, 1> m1{1, 2, 3};
	MatInstance<double, 3, 1> m2{4, 5, 6};
	MatInstance<double, 3, 1> m3{7, 8, 9};
	Armor::Array<double> arr(3, 1, 3);
	arr.get<3>(0) = m1;
	arr.get<3>(1) = m2;
	arr.get<3>(2) = m3;
//	EXPECT_TRUE(m1 == arr.arma_mat(0));
//	EXPECT_TRUE(m1 == arr.arma_vec(0));
//	EXPECT_TRUE(m2 == arr.arma_mat(1));
//	EXPECT_TRUE(m2 == arr.arma_vec(1));
//	EXPECT_TRUE(m3 == arr.arma_mat(2));
//	EXPECT_TRUE(m3 == arr.arma_vec(2));
}
*/


//void fn_armor_mat_fixed(const Armor::Mat<double, 3, 3> &x)
//{}
//
//
//template<class Type, int nr, int nc>
//void fn_armor_mat(const Armor::Mat<Type, nr, nc> &x)
//{}
//
//template<class Type, int nr, int nc>
//void fn_arma_mat(const typename arma::Mat<Type>::template fixed<nr, nc> &x)
//{}
//
//template<class Type, int nr>
//void fn_arma_col(const typename arma::Col<Type>::template fixed<nr> &x)
//{}
//
//template<class Type>
//void fn_arma_scalar(const Type &x)
//{}
//
//
//TEST(Armor_test, conversions) {
//    // Unified Mat representation.
//    MatInstance<double, 3, 3> armor_mat{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//    MatInstance<double, 3, 1> armor_vec{1,2,3};
//    MatInstance<double, 1, 1> armor_scalar{1};
//
//    // Various arma objects
//    arma::Mat<double>::template fixed<3, 3> arma_mat;
//    arma::Mat<double>::template fixed<3, 1> arma_mat_vec;
//    arma::Col<double>::template fixed<3> arma_col_vec;
//    arma::Mat<double>::template fixed<1, 1> arma_mat_scalar;
//    arma::Col<double>::template fixed<1> arma_col_scalar;
//    double scalar;
//
//    // conversions
//    fn_armor_mat_fixed(arma_mat);
//    fn_armor_mat(arma_mat);
//    fn_armor_mat(arma_mat_vec);
//    fn_armor_mat(arma_col_vec);
//    fn_armor_mat(arma_mat_scalar);
//    fn_armor_mat(arma_col_scalar);
//
//
//
//
//}

