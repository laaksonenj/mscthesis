#include <gtest/gtest.h>

#include "fem/math/polynomial/Polynomial1D.hpp"

namespace fem::ut
{
TEST(Polynomial1DTest, Addition)
{
    Polynomial1D p;
    p += Polynomial1D("-t^2+t-1");
    EXPECT_EQ(p, Polynomial1D("-t^2+t-1"));
    p += Polynomial1D("-t");
    EXPECT_EQ(p, Polynomial1D("-t^2-1"));
    p += Polynomial1D("-t^3");
    EXPECT_EQ(p, Polynomial1D("-t^3-t^2-1"));
    p += Monomial1D(mpq_class("21/13"), 14);
    EXPECT_EQ(p, Polynomial1D("21/13t^14-t^3-t^2-1"));
    p = p + 1;
    EXPECT_EQ(p, Polynomial1D("21/13t^14-t^3-t^2"));
    p = p + mpq_class("4/3");
    EXPECT_EQ(p, Polynomial1D("21/13t^14-t^3-t^2+4/3"));
    p += p;
    EXPECT_EQ(p, Polynomial1D("42/13t^14-2t^3-2t^2+8/3"));
}

TEST(Polynomial1DTest, Subtraction)
{
    Polynomial1D p;
    p -= Polynomial1D("t^2-t+1");
    EXPECT_EQ(p, Polynomial1D("-t^2+t-1"));
    p -= Polynomial1D("t");
    EXPECT_EQ(p, Polynomial1D("-t^2-1"));
    p -= Polynomial1D("t^3");
    EXPECT_EQ(p, Polynomial1D("-t^3-t^2-1"));
    p -= Monomial1D(mpq_class("-21/13"), 14);
    EXPECT_EQ(p, Polynomial1D("21/13t^14-t^3-t^2-1"));
    p = p - (-1);
    EXPECT_EQ(p, Polynomial1D("21/13t^14-t^3-t^2"));
    p = p - mpq_class("-4/3");
    EXPECT_EQ(p, Polynomial1D("21/13t^14-t^3-t^2+4/3"));
    p -= p;
    EXPECT_EQ(p, Polynomial1D(0));
}

TEST(Polynomial1DTest, Multiplication)
{
    Polynomial1D p("1");
    p *= Polynomial1D("2-t");
    EXPECT_EQ(p, Polynomial1D("2-t"));
    p = p * p;
    EXPECT_EQ(p, Polynomial1D("4-4t+t^2"));
    p *= -1;
    EXPECT_EQ(p, Polynomial1D("-4+4t-t^2"));
    p = p * Monomial1D(mpq_class("3/2"), 1);
    EXPECT_EQ(p, Polynomial1D("-6t + 6t^2 - 3/2t^3"));
    p *= 0;
    EXPECT_EQ(p, Polynomial1D(0));
}

TEST(Polynomial1DTest, UnaryMinus)
{
    EXPECT_EQ(-Polynomial1D("t^2-1"), Polynomial1D("1-t^2"));
    EXPECT_EQ(-Polynomial1D(4), Polynomial1D(-4));
    EXPECT_EQ(-Polynomial1D(0), Polynomial1D(0));
}

TEST(Polynomial1DTest, ToString)
{
    EXPECT_EQ(toString(Polynomial1D("-t + 3/8t^4 + t + 2t + 0")), "3/8t^4+2t");
    EXPECT_EQ(toString(Polynomial1D("t")), "t");
    EXPECT_EQ(toString(Polynomial1D(0)), "0");
}
} // namespace fem::ut
