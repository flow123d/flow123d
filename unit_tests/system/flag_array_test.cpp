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

    flags_.set(bachelor); // 001
    //std::cout << flags_ << std::endl;
    EXPECT_TRUE( flags_.match(FlagArrayFixture::man_flag) );
    EXPECT_FALSE( flags_.match(FlagArrayFixture::married_flag) );
    EXPECT_FALSE( flags_.match(this->young_flag) );

    flags_.set(this->young_flag);

    Flags::Mask young_bachelor = bachelor & young_flag;
    EXPECT_TRUE(flags_.match(young_bachelor));
    EXPECT_TRUE(flags_.match(bachelor));
    flags_.set(~man_flag);

    EXPECT_FALSE( flags_.match(FlagArrayFixture::man_flag) );
    EXPECT_FALSE( flags_.match(FlagArrayFixture::married_flag) );
    EXPECT_TRUE( flags_.match(this->young_flag) );
    EXPECT_TRUE( flags_.match(this->young_flag & ~man_flag));
    EXPECT_TRUE( flags_.match( young_bachelor & ~man_flag));
}
