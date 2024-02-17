#include <gtest/gtest.h>

#include "fem/math/AffineMap.hpp"

namespace fem::ut
{
TEST(AffineMapTest, Evaluation)
{
    Matrix2mpq A;
    A.col(0) = Vector2mpq("1", "4/3");
    A.col(1) = Vector2mpq("-1", "1");
    Vector2mpq b("-1", "-1");
    AffineMap F(A, b);
    EXPECT_EQ(F({"1", "0"}), Vector2mpq("0", "1/3"));
    EXPECT_EQ(F({"0", "1"}), Vector2mpq("-2", "0"));
    EXPECT_EQ(F({"0", "0"}), Vector2mpq("-1", "-1"));
}

TEST(AffineMapTest, Inverse)
{
    Matrix2mpq A;
    A.col(0) = Vector2mpq("1", "4/3");
    A.col(1) = Vector2mpq("-1", "1");
    Vector2mpq b("-1", "-1");
    AffineMap F(A, b);
    AffineMap Finv = F.inverse();
    EXPECT_EQ(Finv({"0", "1/3"}), Vector2mpq("1", "0"));
    EXPECT_EQ(Finv({"-2", "0"}), Vector2mpq("0", "1"));
    EXPECT_EQ(Finv({"-1", "-1"}), Vector2mpq("0", "0"));
}

TEST(AffineMapTest, Composition)
{
    /* F(x)=Ax+b */
    Matrix2mpq A;
    A.col(0) = Vector2mpq("-1/2", "-1");
    A.col(1) = Vector2mpq("-3/4", "0");
    Vector2mpq b("3/2", "1");
    AffineMap F(A, b);

    /* G(x)=Bx+c */
    Matrix2mpq B;
    B.col(0) = Vector2mpq("1/2", "0");
    B.col(1) = Vector2mpq("0", "1/2");
    Vector2mpq c("1/2", "1/2");
    AffineMap G(B, c);

    /* F(G(x)) */
    Matrix2mpq C;
    C.col(0) = Vector2mpq("-1/4", "-1/2");
    C.col(1) = Vector2mpq("-3/8", "0");
    Vector2mpq d("7/8", "1/2");
    AffineMap FG(C, d);

    EXPECT_EQ(compose(F, G), FG);
}
} // namespace fem::ut
