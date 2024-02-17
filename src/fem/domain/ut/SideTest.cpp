#include <gtest/gtest.h>

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
} // namespace fem::ut
