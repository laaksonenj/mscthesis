#include <gtest/gtest.h>

#include <cmath>

#include "fem/multiprecision/Arithmetic.hpp"

namespace fem::ut
{
TEST(MpArithmeticTest, MpqClassPow)
{
    EXPECT_EQ(pow(mpq_class("-3/4"), -3), mpq_class("-64/27"));
    EXPECT_EQ(pow(mpq_class("-3/4"), -2), mpq_class("16/9"));
    EXPECT_EQ(pow(mpq_class("-3/4"), -1), mpq_class("-4/3"));
    EXPECT_EQ(pow(mpq_class("-3/4"), 0), mpq_class("1"));
    EXPECT_EQ(pow(mpq_class("-3/4"), 1), mpq_class("-3/4"));
    EXPECT_EQ(pow(mpq_class("-3/4"), 2), mpq_class("9/16"));
    EXPECT_EQ(pow(mpq_class("-3/4"), 3), mpq_class("-27/64"));
    EXPECT_EQ(pow(mpq_class("1/2"), 1), mpq_class("1/2"));
    EXPECT_EQ(pow(mpq_class("7/4"), 2), mpq_class("49/16"));
    EXPECT_EQ(pow(mpq_class("2"), 3), mpq_class("8"));
    EXPECT_EQ(pow(mpq_class(0), 0), mpq_class(1));
    EXPECT_EQ(pow(mpq_class(0), 1), mpq_class(0));
}

TEST(MpArithmeticTest, MpzClassLog)
{
    EXPECT_NEAR(log(mpz_class(2)), std::log(2), 1e-7);
    EXPECT_NEAR(log(mpz_class("9223372036854775807")), 43.66827237527655449317, 1e-7);
    EXPECT_NEAR(log(mpz_class("36893488147419103228")), 45.05456673639644511201, 1e-7);
}

TEST(MpArithmeticTest, MpqClassLog)
{
    EXPECT_NEAR(log(mpq_class("1/2")), std::log(0.5), 1e-7);
    EXPECT_NEAR(log(mpq_class("1/2147483647")), -21.487562596892643304518, 1e-7);
    EXPECT_NEAR(log(mpq_class("1/922372036854775807")), -41.3657250477672116069, 1e-7);
    EXPECT_NEAR(log(mpq_class("1/36893488147419103228")), -45.0545667363964451120, 1e-7);
}
} // namespace fem::ut
