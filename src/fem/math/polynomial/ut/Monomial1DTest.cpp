#include <gtest/gtest.h>

#include "fem/math/polynomial/Monomial1D.hpp"

namespace fem::ut
{
TEST(Monomial1DTest, Multiplication)
{
    EXPECT_EQ(Monomial1D(mpq_class("1/4"), 1) * Monomial1D(mpq_class("-4"), 12), Monomial1D(mpq_class(-1), 13));
    EXPECT_EQ(Monomial1D(mpq_class(1), 0) * Monomial1D(mpq_class("-4"), 12), Monomial1D(mpq_class("-4"), 12));
    EXPECT_EQ(Monomial1D(mpq_class(0), 10) * Monomial1D(mpq_class("-4"), 12), Monomial1D(mpq_class(0), 10));
}

TEST(Monomial1DTest, UnaryMinus)
{
    EXPECT_EQ(-Monomial1D(mpq_class("1/4"), 1), Monomial1D(mpq_class("-1/4"), 1));
    EXPECT_EQ(-Monomial1D(mpq_class(0), 1), Monomial1D(mpq_class(0), 1));
}

TEST(Monomial1DTest, ToString)
{
    EXPECT_EQ(toString(Monomial1D(mpq_class("1/4"), 2)), "1/4t^2");
    EXPECT_EQ(toString(Monomial1D(mpq_class("1"), 1)), "t");
    EXPECT_EQ(toString(Monomial1D(mpq_class("-1"), 3)), "-t^3");
    EXPECT_EQ(toString(Monomial1D(mpq_class("0"), 2)), "0");
    EXPECT_EQ(toString(Monomial1D(mpq_class("4"), 0)), "4");
    EXPECT_EQ(toString(Monomial1D(mpq_class("-1/4"), 2024)), "-1/4t^2024");
}
} // namespace fem::ut
