#include <gtest/gtest.h>

#include "fem/math/polynomial/PolynomialComposition.hpp"

namespace fem::ut
{
TEST(PolynomialCompositionTest, Composition1D1D)
{
    EXPECT_EQ(compose(Polynomial1D("t^2"), Polynomial1D("1-t")), Polynomial1D("1-2t+t^2"));
    EXPECT_EQ(compose(Polynomial1D("t+1"), Polynomial1D("t^2-1")), Polynomial1D("t^2"));
    EXPECT_EQ(compose(Polynomial1D("t^2024-t^2-t+1"), Polynomial1D("-t")), Polynomial1D("t^2024-t^2+t+1"));
    EXPECT_EQ(compose(Polynomial1D(4), Polynomial1D("t^3+t")), Polynomial1D(4));
    EXPECT_EQ(compose(Polynomial1D("t^3+t"), Polynomial1D(4)), Polynomial1D(68));
}

TEST(PolynomialCompositionTest, Composition1D2D)
{
    EXPECT_EQ(compose(Polynomial1D("t^2"), Polynomial2D("1-x")), Polynomial2D("1-2x+x^2"));
    EXPECT_EQ(compose(Polynomial1D("t+1"), Polynomial2D("y^2-x-1")), Polynomial2D("y^2-x"));
    EXPECT_EQ(compose(Polynomial1D("t^2024-t^2-t+1"), Polynomial2D("-x")), Polynomial2D("x^2024-x^2+x+1"));
    EXPECT_EQ(compose(Polynomial1D(4), Polynomial2D("x^3+y")), Polynomial2D(4));
    EXPECT_EQ(compose(Polynomial1D("t^3+t"), Polynomial2D(4)), Polynomial2D(68));
}

TEST(PolynomialCompositionTest, Composition2D2D2D)
{
    EXPECT_EQ(compose(Polynomial2D("-3/2x^3y+3/2x^3+3/2xy-3/2x"), Polynomial2D("-x"), Polynomial2D("y")), Polynomial2D("3/2x^3y-3/2x^3-3/2xy+3/2x"));
    EXPECT_EQ(compose(Polynomial2D("xy^2"), Polynomial2D("1-y"), Polynomial2D("x^2")), Polynomial2D("x^4-x^4y"));
    EXPECT_EQ(compose(Polynomial2D("x^2-2xy+y^2"), Polynomial2D(1), Polynomial2D("x")), Polynomial2D("1-2x+x^2"));
    EXPECT_EQ(compose(Polynomial2D("x^2-2xy+y^2"), Polynomial2D(1), Polynomial2D(1)), Polynomial2D(0));
}

TEST(PolynomialCompositionTest, Composition2D1D1D)
{
    EXPECT_EQ(compose(Polynomial2D("-3/2x^3y+3/2x^3+3/2xy-3/2x"), Polynomial1D("t"), Polynomial1D("1-t")), Polynomial1D("3/2t^4 - 3/2t^2"));
    EXPECT_EQ(compose(Polynomial2D("x^2-2xy+y^2"), Polynomial1D(1), Polynomial1D("t")), Polynomial1D("1-2t+t^2"));
    EXPECT_EQ(compose(Polynomial2D("x^2-2xy+y^2"), Polynomial1D(1), Polynomial1D(1)), Polynomial1D(0));
}

TEST(PolynomialCompositionTest, Composition2DAffine)
{
    Matrix2mpq A;
    A << mpq_class("-1/4"), mpq_class("7/2"),
         mpq_class(3), mpq_class(4);
    Vector2mpq b{mpq_class("-1/2"), mpq_class(7)};
    AffineMap G(A, b);
    EXPECT_EQ(compose(Polynomial2D("x^2-2xy+y^2"), G), Polynomial2D("169/16x^2 + 13/4xy + 195/4x + 1/4y^2 + 15/2y + 225/4"));

    A << mpq_class("-1/4"), mpq_class("7/2"),
         mpq_class(3), mpq_class(0);
    b << mpq_class(0), mpq_class(7);
    G = AffineMap(A, b);
    EXPECT_EQ(compose(Polynomial2D("x^2-2xy+y^2"), G), Polynomial2D("169/16x^2 - 91/4xy + 91/2x + 49/4y^2 - 49y + 49"));
}
} // namespace fem::ut
