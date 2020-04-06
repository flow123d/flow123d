#include <flow_gtest.hh>
#include "system/armor.hh"

TEST(Armor_test, indexing1) {
    double mem[4];
	Armor::Mat<double, 3, 1> m1(mem);
	m1[0] = 1;
	m1[1] = 2;
	m1[2] = 3;
	EXPECT_EQ(1, m1[0]);
	EXPECT_EQ(2, m1[1]);
	EXPECT_EQ(3, m1[2]);
}

TEST(Armor_test, equality_operator) {
    double mem[4];
	Armor::Mat<double, 3, 1> m1(mem);
	arma::Mat<double>::fixed<3, 1> a1;
	a1[0] = m1[0] = 1;
	a1[1] = m1[1] = 2;
	a1[2] = m1[2] = 3;
	EXPECT_TRUE(m1 == a1);
	a1[1] = 7;
	EXPECT_FALSE(m1 == a1);
}

//TEST(Armor_test, constructor_list1) {
//	Armor::Mat<double, 3, 1> m1{1, 2, 3};
//	EXPECT_EQ(1, m1[0]);
//	EXPECT_EQ(2, m1[1]);
//	EXPECT_EQ(3, m1[2]);
//}
//
//TEST(Armor_test, constructor_list2) {
//	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
//	EXPECT_TRUE(m1 == a1);
//}
//
//TEST(Armor_test, constructor_arma) {
//	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
//	Armor::Mat<double, 3, 3> m1{a1};
//	EXPECT_TRUE(m1 == a1);
//}

TEST(Armor_test, indexing2) {
	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
	double mem[9];
	Armor::Mat<double, 3, 3> m1(mem);
	for (uint i{0}; i < 3; ++i) {
		for (uint j{0}; j < 3; ++j) {
			m1(i,j) = a1(i,j);
		}
	}
	EXPECT_TRUE(m1 == a1);
}

TEST(Armor_test, assignment_list1) {
    double mem[3];
	Armor::Mat<double, 3, 1> m1(mem);
	m1 = {1, 2, 3};
	EXPECT_EQ(1, m1[0]);
	EXPECT_EQ(2, m1[1]);
	EXPECT_EQ(3, m1[2]);
}

TEST(Armor_test, assignment_list2) {
	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
	double mem[9];
	Armor::Mat<double, 3, 3> m1(mem);
	m1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	EXPECT_TRUE(m1 == a1);
}

TEST(Armor_test, assignment_armor) {
    double mem1[9], mem2[9];
	Armor::Mat<double, 3, 3> m1(mem1), m2(mem2);
	m2 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	m1 = m2;
	EXPECT_TRUE(m1 == m2);
}

TEST(Armor_test, assignment_arma) {
    double mem[9];
	Armor::Mat<double, 3, 3> m1(mem);
	arma::Mat<double>::fixed<3, 3> a1{1, 4, 7, 2, 5, 8, 3, 6, 9};
	m1 = a1;
	EXPECT_TRUE(m1 == a1);
}

TEST(Armor_test, size) {
    double mem1[3], mem2[9];
	Armor::Mat<double, 3, 1> m1(mem1);
	Armor::Mat<double, 3, 3> m2(mem2);
	EXPECT_EQ(3, m1.size());
	EXPECT_EQ(9, m2.size());
}

TEST(Armor_test, memptr) {
    double mem1[9];
    Armor::Mat<double, 3, 3> m1(mem1);
    m1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
	EXPECT_TRUE(m1 == a1);
}

//TEST(Armor_test, arma) {
//	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
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
//    Armor::Mat<double, 3, 3> m1(mem1), m2(mem2);
//    m2 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//
//	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
//	double a = dot(a1, a1);
//	double m = dot(m1, m1);
//	EXPECT_EQ(a, m);
//}
//
//TEST(Armor_test, plus) {
//	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
//	a1 = a1 + a1;
//	m1 = m1 + m1;
//	EXPECT_TRUE(m1 == a1);
//}
//
//TEST(Armor_test, minus) {
//	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
//	a1 = a1 - a1;
//	m1 = m1 - m1;
//	EXPECT_TRUE(m1 == a1);
//}
//
//TEST(Armor_test, multiplication) {
//	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
//	a1 = a1 * a1;
//	m1 = m1 * m1;
//	EXPECT_TRUE(m1 == a1);
//}
//
//TEST(Armor_test, multiplication_per_elements) {
//	Armor::Mat<double, 3, 3> m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//	arma::Mat<double>::fixed<3, 3> a1{m1.memptr()};
//	a1 = a1 % a1;
//	m1 = m1 % m1;
//	EXPECT_TRUE(m1 == a1);
//}
