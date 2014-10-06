/*
 * flag_array.cpp
 *
 *  Created on: May 2, 2014
 *      Author: jb
 */

#include <flow_gtest.hh>
#include <system/flag_array.hh>
#include <iostream>

class FlagArrayFixture : public testing::Test
{
public:
    FlagArrayFixture()  {

    }
    ~FlagArrayFixture() {}

    typedef FlagArray<FlagArrayFixture, 4> Flags;
    Flags  flags_;

    static constexpr Flags::Mask man_flag   {1};
    static constexpr Flags::Mask married_flag   {2};
    static constexpr Flags::Mask bachelor   {man_flag & ~married_flag};


    const Flags::Mask young_flag {4};
    const Flags::Mask old_maid   {~man_flag & ~married_flag & ~young_flag};
    const Flags::Mask maid   {~man_flag & ~married_flag};

};

constexpr FlagArrayFixture::Flags::Mask FlagArrayFixture::man_flag;
constexpr FlagArrayFixture::Flags::Mask FlagArrayFixture::married_flag;
constexpr FlagArrayFixture::Flags::Mask FlagArrayFixture::bachelor;

TEST_F(FlagArrayFixture, set_and_test) {

    EXPECT_FALSE( flags_.match(FlagArrayFixture::man_flag) );
    EXPECT_FALSE( flags_.match(FlagArrayFixture::married_flag) );
    EXPECT_FALSE( flags_.match(this->young_flag) );

    //std::cout << man_flag << std::endl;
    //std::cout << ~man_flag << std::endl;
    //std::cout << bachelor << std::endl;

    flags_.add(bachelor); // 001
    //std::cout << flags_ << std::endl;
    EXPECT_TRUE( flags_.match(FlagArrayFixture::man_flag) );
    EXPECT_FALSE( flags_.match(FlagArrayFixture::married_flag) );
    EXPECT_FALSE( flags_.match(this->young_flag) );

    flags_.add(this->young_flag);

    Flags::Mask young_bachelor = bachelor & young_flag;
    EXPECT_TRUE(flags_.match(young_bachelor));
    EXPECT_TRUE(flags_.match(bachelor));
    flags_.add(~man_flag);

    EXPECT_FALSE( flags_.match(FlagArrayFixture::man_flag) );
    EXPECT_FALSE( flags_.match(FlagArrayFixture::married_flag) );
    EXPECT_TRUE( flags_.match(this->young_flag) );
    EXPECT_TRUE( flags_.match(this->young_flag & ~man_flag));
    EXPECT_TRUE( flags_.match( young_bachelor & ~man_flag));

    // operator ==
    bool result= ( Flags(this->young_flag)==flags_ );
}

void tst(int f1, int f2, int m1, int m2, int s1, int s2, int r) {
    typedef FlagArray<FlagArrayFixture, 4> Flags;
    Flags::Mask m(2*m2+m1, 2*s2+s1);
    Flags f(2*f2+f1);

    //std::cout << f1<<f2<<m1<<m2<<s1<<s2<<std::endl;
    if (r) EXPECT_TRUE(f.match(m));
    else EXPECT_FALSE(f.match(m));
}
TEST_F(FlagArrayFixture, match) {

    tst(0,0,0,0,0,0,1);
    tst(0,0,0,0,0,1,1);
    tst(0,0,0,0,1,0,1);
    tst(0,0,0,0,1,1,1);
    tst(0,0,0,1,0,0,1);
    tst(0,0,0,1,0,1,0);
    tst(0,0,0,1,1,0,1);
    tst(0,0,0,1,1,1,0);
    tst(0,0,1,0,0,0,1);
    tst(0,0,1,0,0,1,1);
    tst(0,0,1,0,1,0,0);
    tst(0,0,1,0,1,1,0);
    tst(0,0,1,1,0,0,1);
    tst(0,0,1,1,0,1,0);
    tst(0,0,1,1,1,0,0);
    tst(0,0,1,1,1,1,0);

    tst(0,1,0,0,0,0,1);
    tst(0,1,0,0,0,1,1);
    tst(0,1,0,0,1,0,1);
    tst(0,1,0,0,1,1,1);
    tst(0,1,0,1,0,0,0);
    tst(0,1,0,1,0,1,1);
    tst(0,1,0,1,1,0,0);
    tst(0,1,0,1,1,1,1);
    tst(0,1,1,0,0,0,1);
    tst(0,1,1,0,0,1,1);
    tst(0,1,1,0,1,0,0);
    tst(0,1,1,0,1,1,0);
    tst(0,1,1,1,0,0,0);
    tst(0,1,1,1,0,1,1);
    tst(0,1,1,1,1,0,0);
    tst(0,1,1,1,1,1,0);

    tst(1,0,0,0,0,0,1);
    tst(1,0,0,0,0,1,1);
    tst(1,0,0,0,1,0,1);
    tst(1,0,0,0,1,1,1);
    tst(1,0,0,1,0,0,1);
    tst(1,0,0,1,0,1,0);
    tst(1,0,0,1,1,0,1);
    tst(1,0,0,1,1,1,0);
    tst(1,0,1,0,0,0,0);
    tst(1,0,1,0,0,1,0);
    tst(1,0,1,0,1,0,1);
    tst(1,0,1,0,1,1,1);
    tst(1,0,1,1,0,0,0);
    tst(1,0,1,1,0,1,0);
    tst(1,0,1,1,1,0,1);
    tst(1,0,1,1,1,1,0);

    tst(1,1,0,0,0,0,1);
    tst(1,1,0,0,0,1,1);
    tst(1,1,0,0,1,0,1);
    tst(1,1,0,0,1,1,1);
    tst(1,1,0,1,0,0,0);
    tst(1,1,0,1,0,1,1);
    tst(1,1,0,1,1,0,0);
    tst(1,1,0,1,1,1,1);
    tst(1,1,1,0,0,0,0);
    tst(1,1,1,0,0,1,0);
    tst(1,1,1,0,1,0,1);
    tst(1,1,1,0,1,1,1);
    tst(1,1,1,1,0,0,0);
    tst(1,1,1,1,0,1,0);
    tst(1,1,1,1,1,0,0);
    tst(1,1,1,1,1,1,1);
}
