#include <gtest/gtest.h>

#include "fem/math/polynomial/Monomial2D.hpp"

namespace fem::ut
{
TEST(Monomial2DTest, Multiplication)
{
    EXPECT_EQ(Monomial2D(mpq_class("-5/7"), 1, 4) * Monomial2D(mpq_class("1/2"), 10, 3), Monomial2D(mpq_class("-5/14"), 11, 7));
    EXPECT_EQ(Monomial2D(mpq_class(1), 0, 4) * Monomial2D(mpq_class("-3/8"), 10, 0), Monomial2D(mpq_class("-3/8"), 10, 4));
    EXPECT_EQ(Monomial2D(mpq_class(0), 0, 4) * Monomial2D(mpq_class("-3/8"), 10, 0), Monomial2D(mpq_class(0), 10, 4));
}

TEST(Monomial2DTest, UnaryMinus)
{
    EXPECT_EQ(-Monomial2D(mpq_class("-5/7"), 1, 4), Monomial2D(mpq_class("5/7"), 1, 4));
    EXPECT_EQ(-Monomial2D(mpq_class("5/7"), 1, 4), Monomial2D(mpq_class("-5/7"), 1, 4));
    EXPECT_EQ(-Monomial2D(mpq_class(0), 0, 4), Monomial2D(mpq_class(0), 0, 0));
}

TEST(Monomial2DTest, Evaluation)
{
    EXPECT_EQ(Monomial2D(mpq_class(1), 2, 2)({mpq_class("1/2"), mpq_class("1/3")}), mpq_class("1/36"));
    EXPECT_EQ(Monomial2D(mpq_class("5/4"), 0, 3)({mpq_class("1/4"), mpq_class(2)}), mpq_class(10));
    EXPECT_EQ(Monomial2D(mpq_class("5/4"), 0, 3)({mpq_class(0), mpq_class(2)}), mpq_class(10));
    EXPECT_EQ(Monomial2D(mpq_class(0), 1, 3)({mpq_class("1/4"), mpq_class(2)}), mpq_class(0));
    EXPECT_EQ(Monomial2D(mpq_class(4), 1, 3)({mpq_class("1/4"), mpq_class(0)}), mpq_class(0));
}

TEST(Monomial2DTest, ToString)
{
    EXPECT_EQ(toString(Monomial2D(mpq_class("1/4"), 2, 3)), "1/4x^2y^3");
    EXPECT_EQ(toString(Monomial2D(mpq_class("1"), 1, 3)), "xy^3");
    EXPECT_EQ(toString(Monomial2D(mpq_class("-1"), 3, 2024)), "-x^3y^2024");
    EXPECT_EQ(toString(Monomial2D(mpq_class("0"), 2, 2)), "0");
    EXPECT_EQ(toString(Monomial2D(mpq_class("4"), 0, 0)), "4");
    EXPECT_EQ(toString(Monomial2D(mpq_class("4"), 1, 0)), "4x");
}
} // namespace fem::ut
