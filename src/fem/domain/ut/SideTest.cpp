#include <gtest/gtest.h>

#include <cmath>

#include "fem/domain/Side.hpp"

namespace fem::ut
{
TEST(SideTest, Equality)
{
    const Side side1(Node(mpq_class(0), mpq_class(0)), Node(mpq_class("-1/7"), mpq_class("7/3")));
    const Side side2(Node(mpq_class("-1/7"), mpq_class("7/3")), Node(mpq_class(0), mpq_class(0)));
    const Side side3(Node(mpq_class(0), mpq_class(-1)), Node(mpq_class(0), mpq_class(1)));
    EXPECT_EQ(side1, side1);
    EXPECT_EQ(side1, side2);
    EXPECT_EQ(side2, side1);
    EXPECT_NE(side1, side3);
    EXPECT_NE(side3, side1);
}

TEST(SideTest, AreParallel)
{
    const Side side1(Node(mpq_class(0), mpq_class(0)), Node(mpq_class("-1/7"), mpq_class("7/3")));
    const Side side2(Node(mpq_class("-1/7"), mpq_class("7/3")), Node(mpq_class(0), mpq_class(0)));
    const Side side3(Node(mpq_class(0), mpq_class(-1)), Node(mpq_class(0), mpq_class(1)));
    const Side side4(Node(mpq_class(7), mpq_class(128)), Node(mpq_class(7), mpq_class("78/6")));
    EXPECT_TRUE(areParallel(side1, side1));
    EXPECT_TRUE(areParallel(side1, side2));
    EXPECT_TRUE(areParallel(side2, side1));
    EXPECT_TRUE(areParallel(side3, side4));
    EXPECT_FALSE(areParallel(side1, side3));
    EXPECT_FALSE(areParallel(side3, side1));
    EXPECT_FALSE(areParallel(side2, side4));
}

TEST(SideTest, Intersection)
{
    const Side side1(Node(mpq_class(0), mpq_class(0)), Node(mpq_class("-1/7"), mpq_class("7/3")));
    const Side side2(Node(mpq_class("-1/7"), mpq_class("7/3")), Node(mpq_class(0), mpq_class(0)));
    const Side side3(Node(mpq_class(0), mpq_class(-1)), Node(mpq_class(0), mpq_class(1)));
    const Side side4(Node(mpq_class(7), mpq_class(128)), Node(mpq_class(7), mpq_class("78/6")));
    const Side side5(Node(mpq_class(0), mpq_class(1)), Node(mpq_class(10), mpq_class(50)));
    EXPECT_TRUE(areIntersecting(side1, side1));
    EXPECT_TRUE(areIntersecting(side1, side2));
    EXPECT_TRUE(areIntersecting(side1, side3));
    EXPECT_TRUE(areIntersecting(side3, side1));
    EXPECT_TRUE(areIntersecting(side3, side2));
    EXPECT_TRUE(areIntersecting(side3, side5));
    EXPECT_TRUE(areIntersecting(side4, side5));
    EXPECT_FALSE(areIntersecting(side1, side4));
    EXPECT_FALSE(areIntersecting(side2, side4));
    EXPECT_FALSE(areIntersecting(side1, side5));
    EXPECT_FALSE(areIntersecting(side2, side5));
    EXPECT_FALSE(areIntersecting(side3, side4));
}

namespace
{
void checkNear(const Vector2mpq& actual, const Vector2mpq& expected)
{
    EXPECT_NEAR(actual(0).get_d(), expected(0).get_d(), 1e-7);
    EXPECT_NEAR(actual(1).get_d(), expected(1).get_d(), 1e-7);
}
} // namespace

TEST(SideTest, Normal)
{
    checkNear(calculateNormal(Side({1, -1}, {1, 1})), Vector2mpq(1, 0));
    checkNear(calculateNormal(Side({0, 0}, {"-1/4", "1/4"})), Vector2mpq(1.0/std::sqrt(2.0), 1.0/std::sqrt(2.0)));
    checkNear(calculateNormal(Side({"3/8", "12/5"}, {"-4", "0"})), Vector2mpq(-96.0/std::sqrt(39841), 175.0/std::sqrt(39841)));
}
} // namespace fem::ut
