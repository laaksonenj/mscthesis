#include <gtest/gtest.h>

#include "fem/math/polynomial/Polynomial2D.hpp"

namespace fem::ut
{
TEST(Polynomial2DTest, Addition)
{
    Polynomial2D p;
    p += Polynomial2D("x^2+2xy+y^2");
    EXPECT_EQ(p, Polynomial2D("x^2+2xy+y^2"));
    p += Monomial2D(mpq_class("-2"), 1, 1);
    EXPECT_EQ(p, Polynomial2D("x^2+y^2"));
    p = p + mpq_class("3/7");
    EXPECT_EQ(p, Polynomial2D("x^2+y^2+3/7"));
    p = p + 1;
    EXPECT_EQ(p, Polynomial2D("x^2+y^2+10/7"));
    p += p;
    EXPECT_EQ(p, Polynomial2D("2x^2+2y^2+20/7"));
}

TEST(Polynomial2DTest, Subtraction)
{
    Polynomial2D p;
    p -= Polynomial2D("2xy - 4/7");
    EXPECT_EQ(p, Polynomial2D("-2xy+4/7"));
    p = p - Monomial2D(mpq_class(1), 0, 2);
    EXPECT_EQ(p, Polynomial2D("-y^2-2xy+4/7"));
    p -= mpq_class("4/7");
    EXPECT_EQ(p, Polynomial2D("-y^2-2xy"));
    p -= 1;
    EXPECT_EQ(p, Polynomial2D("-y^2-2xy-1"));
    p = p - p;
    EXPECT_EQ(p, Polynomial2D(0));
}

TEST(Polynomial2DTest, Multiplication)
{
    Polynomial2D p = 1;
    p *= Polynomial2D("x^2+2xy+y^2");
    EXPECT_EQ(p, Polynomial2D("x^2+2xy+y^2"));
    p = p * Polynomial2D("4/7x^7 - 8/3");
    EXPECT_EQ(p, Polynomial2D("4/7x^9-8/3x^2+8/7x^8y-16/3xy+4/7x^7y^2-8/3y^2"));
    EXPECT_EQ(1 * p, p);
    EXPECT_EQ(0 * p, 0);
}

TEST(Polynomial2DTest, UnaryMinus)
{
    EXPECT_EQ(-Polynomial2D("x^2-2xy+y^2"), Polynomial2D("-x^2+2xy-y^2"));
    EXPECT_EQ(-Polynomial2D(4), Polynomial2D(-4));
    EXPECT_EQ(-Polynomial2D(0), Monomial2D(0, 1, 1));
}

TEST(Polynomial2DTest, Evaluation)
{
    EXPECT_EQ(Polynomial2D("x^2y^2")({mpq_class("1/2"), mpq_class("1/3")}), mpq_class("1/36"));
    EXPECT_EQ(Polynomial2D("-y^7 x^3 + 1 + x")({mpq_class("3/7"), mpq_class("1")}), mpq_class("463/343"));
    EXPECT_EQ(Polynomial2D("x^11- x")({mpq_class("1/2"), mpq_class("123/456")}), mpq_class("-1023/2048"));
    EXPECT_EQ(Polynomial2D("x^2+2xy+y^2")({mpq_class("-4/3"), mpq_class("0")}), mpq_class("16/9"));
    EXPECT_EQ(Polynomial2D("  4 + 0 + x")({mpq_class(-4), mpq_class("1/5")}), mpq_class(0));
}
} // namespace fem::ut
