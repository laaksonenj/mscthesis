#include <gtest/gtest.h>

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
    EXPECT_EQ(pow(mpq_class(0), 0), mpq_class(1));
    EXPECT_EQ(pow(mpq_class(0), 1), mpq_class(0));
}
} // namespace fem::ut
