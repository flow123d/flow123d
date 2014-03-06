#include <flow_gtest.hh>

#include "mesh/bounding_box.hh"


TEST(BoundingBox, constructors) {
	typedef BoundingBox::Point Point;
	BoundingBox box; 	//default

	BoundingBox box1( Point(" 1 2 3"), Point("3 4 5") );  // basic constructor
	EXPECT_TRUE( arma::min( box1.min() == Point("1 2 3") ) );
	EXPECT_TRUE( arma::min( box1.max() == Point("3 4 5") ) );

	BoundingBox box2( box1); // Copy
    EXPECT_TRUE( arma::min( box2.min() == Point("1 2 3") ) );
	EXPECT_TRUE( arma::min( box2.max() == Point("3 4 5") ) );

	vector<Point> b_list = {
			Point("3 3 3"),
			Point("1 4 5"),
			Point("2 -10 4")
	};

	BoundingBox box3(b_list);
	EXPECT_TRUE( arma::min( box3.min() == Point("1 -10 3") ) );
	EXPECT_TRUE( arma::min( box3.max() == Point("3 4 5") ) );

}

class BB :  public testing::Test
{
public:
	BB() :
		box1( Point("1 2 3"), Point("3 2 3")),
		box2( Point("1 2 3"), Point("3 4 3")),
		box3( Point(" 1 2 3"), Point("3 4 5") )
	{}


	~BB() {
	}

	typedef BoundingBox::Point Point;

	BoundingBox box1, box2, box3;

};



TEST_F(BB, access) {

	EXPECT_TRUE( arma::min( box3.center() == Point("2 3 4") ) );
	EXPECT_DOUBLE_EQ( 2, box3.projection_center(0) );
	EXPECT_DOUBLE_EQ( 3, box3.projection_center(1) );
	EXPECT_DOUBLE_EQ( 4, box3.projection_center(2) );

	EXPECT_ASSERT_DEATH( {box3.projection_center(4);}, "");
}

TEST_F(BB, contains_point) {

    EXPECT_TRUE( box3.contains_point(Point("2 3 4")) );
    EXPECT_TRUE( box3.contains_point(Point("1 2 3")) );
    EXPECT_TRUE( box3.contains_point(Point("3 4 5")) );
    EXPECT_FALSE( box3.contains_point(Point("0.999 2 3")) );
    EXPECT_FALSE( box3.contains_point(Point("3 4 5.0001")) );

    EXPECT_TRUE( box2.contains_point(Point("2 3 3")) );
    EXPECT_TRUE( box2.contains_point(Point("1 4 3")) );

    EXPECT_TRUE( box1.contains_point(Point("2 2 3")) );
    EXPECT_TRUE( box1.contains_point(Point("1 2 3")) );
}

TEST_F(BB, intersect) {

    EXPECT_TRUE( box3.intersect( box3) );
    EXPECT_TRUE( box3.intersect( box2) );
    EXPECT_TRUE( box3.intersect( box1) );
    EXPECT_TRUE( box3.intersect(BoundingBox(
    		Point("1.5 2.5 3.5"), Point("3.5 4.5 5.5") ) ) );

    EXPECT_TRUE( box2.intersect( box3) );
    EXPECT_TRUE( box2.intersect( box2) );
    EXPECT_TRUE( box2.intersect( box1) );
    EXPECT_TRUE( box2.intersect(BoundingBox(
    		Point("1.5 2.5 3"), Point("3.5 4.5 3") ) ) );

    EXPECT_TRUE( box1.intersect( box3) );
    EXPECT_TRUE( box1.intersect( box2) );
    EXPECT_TRUE( box1.intersect( box1) );
    EXPECT_TRUE( box1.intersect(BoundingBox(
    		Point("1.5 2 3"), Point("3.5 2 3") ) ) );

    EXPECT_FALSE( box3.intersect(BoundingBox(
    		Point("3.5 2 3"), Point("4 4 5") ) ) );
    EXPECT_FALSE( box2.intersect(BoundingBox(
     		Point("3.5 2 0"), Point("4 4 0") ) ) );
    EXPECT_FALSE( box1.intersect(BoundingBox(
     		Point("3.5 0 0"), Point("4 0 0") ) ) );
}

TEST_F(BB, split) {
	BoundingBox b1, b2;

	box3.split(2, 4.5, b1, b2);
	cout << b1;
	cout << b2;
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b1.min() - Point("1 2 3") , "inf") );
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b1.max() - Point("3 4 4.5") , "inf") );
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b2.min() - Point("1 2 4.5") , "inf") );
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b2.max() - Point("3 4 5") , "inf") );

	box2.split(1, 2.5, b1, b2);
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b1.min() - Point("1 2 3") , "inf") );
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b1.max() - Point("3 2.5 3") , "inf") );
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b2.min() - Point("1 2.5 3") , "inf") );
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b2.max() - Point("3 4 3") , "inf") );

	box1.split(0, 2.0, b1, b2);
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b1.min() - Point("1 2 3") , "inf") );
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b1.max() - Point("2 2 3") , "inf") );
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b2.min() - Point("2 2 3") , "inf") );
	EXPECT_DOUBLE_EQ( 0.0, arma::norm( b2.max() - Point("3 2 3") , "inf") );


}

// BoundingBox::expand is tested  as part of constructor from vector of points

TEST_F(BB, longest_axis) {
	box3.expand(Point("3 5 7"));  // (3 5 7) - (1 2 3 ) = (2 3 4)
	EXPECT_EQ( 2 , box3.longest_axis() );
	box3.expand(Point("3 7 5"));  // (3 7 5) - (1 2 3 ) = (2 5 2)
	EXPECT_EQ( 1 , box3.longest_axis() );
	box3.expand(Point("7 5 3"));  // (7 5 3) - (1 2 3 ) = (6 3 0)
	EXPECT_EQ( 0 , box3.longest_axis() );

}
