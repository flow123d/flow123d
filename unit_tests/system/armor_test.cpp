#include <flow_gtest.hh>
#include "system/armor.hh"

TEST(Armor_test, indexing1) {
	armor::Mat<double, 3, 1> m1;
	m1[0] = 1;
	m1[1] = 2;
	m1[2] = 3;
	EXPECT_EQ(1, m1[0]);
	EXPECT_EQ(2, m1[1]);
	EXPECT_EQ(3, m1[2]);
}

TEST(Armor_test, equality_operator) {
	armor::Mat<double, 3, 1> m1;
	arma::Mat<double>::fixed<3, 1> a1;
	a1[0] = m1[0] = 1;
	a1[1] = m1[1] = 2;
	a1[2] = m1[2] = 3;
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
	a1[1] = 7;
	EXPECT_FALSE(approx_equal(m1, a1, "absdiff", 0.0));
}

TEST(Armor_test, constructor_list1) {
	armor::Mat<double, 3, 1> m1{1, 2, 3};
	EXPECT_EQ(1, m1[0]);
	EXPECT_EQ(2, m1[1]);
	EXPECT_EQ(3, m1[2]);
}

TEST(Armor_test, constructor_list2) {
	armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
}

TEST(Armor_test, constructor_arma) {
	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
	armor::Mat<double, 3, 3> m1{a1};
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
}

TEST(Armor_test, indexing2) {
	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
	armor::Mat<double, 3, 3> m1;
	for (uint i{0}; i < 3; ++i) {
		for (uint j{0}; j < 3; ++j) {
			m1(i,j) = a1(i,j);
		}
	}
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
}

TEST(Armor_test, assignment_list1) {
	armor::Mat<double, 3, 1> m1;
	m1 = {1, 2, 3};
	EXPECT_EQ(1, m1[0]);
	EXPECT_EQ(2, m1[1]);
	EXPECT_EQ(3, m1[2]);
}

TEST(Armor_test, assignment_list2) {
	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
	armor::Mat<double, 3, 3> m1;
	m1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
}

TEST(Armor_test, assignment_armor) {
	armor::Mat<double, 3, 3> m1, m2;
	m2 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	m1 = m2;
	EXPECT_TRUE(approx_equal(m1, m2, "absdiff", 0.0));
}

TEST(Armor_test, assignment_arma) {
	armor::Mat<double, 3, 3> m1;
	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
	m1 = a1;
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
}

TEST(Armor_test, size) {
	armor::Mat<double, 3, 1> m1;
	armor::Mat<double, 3, 3> m2;
	EXPECT_EQ((uint) 3, m1.size());
	EXPECT_EQ((uint) 9, m2.size());
}

TEST(Armor_test, memptr) {
	armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
}

TEST(Armor_test, arma) {
	armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
	arma::Mat<double>::fixed<3, 3> a2 = m1.arma();
	EXPECT_TRUE(approx_equal(a1, a2, "absdiff", 0.0));
}

TEST(Armor_test, dot) {
	armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
	double a = dot(a1, a1);
	double m = dot(m1, m1);
	EXPECT_EQ(a, m);
}

TEST(Armor_test, plus) {
	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
	a1 = a1 + a1;
	m1 = m1 + m1;
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
}

TEST(Armor_test, minus) {
	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
	a1 = a1 - a1;
	m1 = m1 - m1;
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
}

TEST(Armor_test, multiplication) {
	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
	a1 = a1 * a1;
	m1 = m1 * m1;
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
}

TEST(Armor_test, multiplication_per_elements) {
	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
	a1 = a1 % a1;
	m1 = m1 % m1;
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
}

TEST(Armor_test, zeros) {
	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
    m1.zeros();
    a1.zeros();
	EXPECT_TRUE(approx_equal(m1, a1, "absdiff", 0.0));
}
