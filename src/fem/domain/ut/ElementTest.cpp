#include <gtest/gtest.h>

#include "fem/domain/Element.hpp"

namespace fem::ut
{
TEST(ElementTest, ParallelogramIsCorrect)
{
    const std::vector<Node> nodes = {Node(mpq_class("3/2"), mpq_class(1)),
                                     Node(mpq_class("3/4"), mpq_class(1)),
                                     Node(mpq_class("1/4"), mpq_class(0)),
                                     Node(mpq_class(1), mpq_class(0))};
    const std::vector<Side> sides = {Side(nodes[0], nodes[1]),
                                     Side(nodes[1], nodes[2]),
                                     Side(nodes[2], nodes[3]),
                                     Side(nodes[3], nodes[0])};
    Matrix2mpq A;
    A.col(0) = Vector2mpq("-3/8", "0");
    A.col(1) = Vector2mpq("-1/4", "-1/2");
    const Vector2mpq b("7/8", "1/2");
    const AffineMap F(A, b);

    const Parallelogram quad(nodes[0], nodes[1], nodes[2], nodes[3]);
    EXPECT_EQ(quad.getElementType(), ElementType_Parallelogram);
    EXPECT_EQ(quad.getNumOfNodes(), 4);
    EXPECT_EQ(quad.getNumOfSides(), 4);
    for (int i = 0; i < 4; i++)
    {
        EXPECT_EQ(quad.getNode(i), nodes[i]);
        EXPECT_EQ(quad.getSide(i), sides[i]);
    }
    EXPECT_EQ(quad.getReferenceElementMap(), F);
}

TEST(ElementTest, TriangleIsCorrect)
{
    const std::vector<Node> nodes = {Node(mpq_class("7/2"), mpq_class("8/3")),
                                     Node(mpq_class(-1), mpq_class(0)),
                                     Node(mpq_class(2), mpq_class(1))};
    const std::vector<Side> sides = {Side(nodes[0], nodes[1]),
                                     Side(nodes[1], nodes[2]),
                                     Side(nodes[2], nodes[0])};
    Matrix2mpq A;
    A.col(0) = Vector2mpq("-9/2", "-8/3");
    A.col(1) = Vector2mpq("-3/2", "-5/3");
    const Vector2mpq b("7/2", "8/3");
    const AffineMap F(A, b);

    const Triangle tri(nodes[0], nodes[1], nodes[2]);
    EXPECT_EQ(tri.getElementType(), ElementType_Triangle);
    EXPECT_EQ(tri.getNumOfNodes(), 3);
    EXPECT_EQ(tri.getNumOfSides(), 3);
    for (int i = 0; i < 3; i++)
    {
        EXPECT_EQ(tri.getNode(i), nodes[i]);
        EXPECT_EQ(tri.getSide(i), sides[i]);
    }
    EXPECT_EQ(tri.getReferenceElementMap(), F);
}

TEST(ElementTest, DegenerateElementAborts)
{
    EXPECT_DEATH(Parallelogram(Node(0, 0), Node(1, 0), Node(2, 0), Node(3, 0)), ".*");
    EXPECT_DEATH(Triangle(Node(0, 0), Node(1, 0), Node(2, 0)), ".*");
}

TEST(ElementTest, ParallelogramAbortsIfSidesAreNotParallel)
{
    EXPECT_DEATH(Parallelogram(Node(0, 0), Node(1, -1), Node(1, 1), Node(0, 1)), ".*");
    EXPECT_DEATH(Parallelogram(Node(0, 0), Node(2, 0), Node(1, 1), Node(0, 1)), ".*");
}

TEST(ElementTest, OrientationMustBeCounterClockwise)
{
    EXPECT_DEATH(Parallelogram(Node(mpq_class("3/2"), 1), Node(1, 0), Node(mpq_class("1/4"), 0), Node(mpq_class("3/4"), 1)), ".*");
    EXPECT_DEATH(Triangle(Node(mpq_class("3/2"), 1), Node(1, 0), Node(mpq_class("1/4"), 0)), ".*");
}

TEST(ElementTest, Intersection)
{
    const Parallelogram el1(Node(0, 0), Node(0, -2), Node(2, -3), Node(2, -1));
    const Parallelogram el2(Node(1, -1), Node(1, -3), Node(2, -3), Node(2, -1));
    const Parallelogram el3(Node(2, -1), Node(2, -3), Node(3, -3), Node(3, -1));
    const Parallelogram el4(Node(4, -1), Node(5, -1), Node(5, 0), Node(4, 0));
    const Parallelogram el5(Node(0, -1), Node(-1, 0), Node(-2, -1), Node(-1, -2));
    const Triangle el6(Node(2, -3), Node(2, -4), Node(3, -4));
    const Triangle el7(Node(2, -2), Node(2, -3), Node(3, -3));
    const Triangle el8(Node(mpq_class("3/2"), -2), Node(1, -1), Node(mpq_class("1/2"), -2));
    const Triangle el9(Node(1, -1), Node(0, 0), Node(mpq_class("1/2"), -1));

    EXPECT_TRUE(areIntersecting(el1, el2));
    EXPECT_TRUE(areIntersecting(el1, el3));
    EXPECT_FALSE(areIntersecting(el1, el4));
    EXPECT_TRUE(areIntersecting(el1, el5));
    EXPECT_TRUE(areIntersecting(el1, el6));
    EXPECT_TRUE(areIntersecting(el1, el7));
    EXPECT_TRUE(areIntersecting(el1, el8));
    EXPECT_TRUE(areIntersecting(el1, el9));

    EXPECT_FALSE(isIntersectionOneNode(el1, el2));
    EXPECT_FALSE(isIntersectionOneNode(el1, el3));
    EXPECT_FALSE(isIntersectionOneNode(el1, el4));
    EXPECT_FALSE(isIntersectionOneNode(el1, el5));
    EXPECT_TRUE(isIntersectionOneNode(el1, el6));
    EXPECT_FALSE(isIntersectionOneNode(el1, el7));
    EXPECT_FALSE(isIntersectionOneNode(el1, el8));
    EXPECT_FALSE(isIntersectionOneNode(el1, el9));

    EXPECT_FALSE(isIntersectionOneSide(el1, el2));
    EXPECT_TRUE(isIntersectionOneSide(el1, el3));
    EXPECT_FALSE(isIntersectionOneSide(el1, el4));
    EXPECT_FALSE(isIntersectionOneSide(el1, el5));
    EXPECT_FALSE(isIntersectionOneSide(el1, el6));
    EXPECT_FALSE(isIntersectionOneSide(el1, el7));
    EXPECT_FALSE(isIntersectionOneSide(el1, el8));
    EXPECT_FALSE(isIntersectionOneSide(el1, el9));
}

TEST(ElementTest, IsPointInsideElement)
{
    const Parallelogram quad(Node(-1, -1), Node(1, -1), Node(1, 1), Node(-1, 1));
    const Triangle tri(Node(0, 0), Node(1, 0), Node(0, 1));
    EXPECT_TRUE(isPointInsideElement({0, 0}, quad));
    EXPECT_TRUE(isPointInsideElement({0, 0}, tri));
    EXPECT_TRUE(isPointInsideElement({-1, 0}, quad));
    EXPECT_FALSE(isPointInsideElement({-1, 0}, tri));
    EXPECT_TRUE(isPointInsideElement({mpq_class("1/4"), mpq_class("1/4")}, quad));
    EXPECT_TRUE(isPointInsideElement({mpq_class("1/4"), mpq_class("1/4")}, tri));
    EXPECT_FALSE(isPointInsideElement({mpq_class("9/4"), mpq_class("9/4")}, quad));
    EXPECT_FALSE(isPointInsideElement({mpq_class("9/4"), mpq_class("9/4")}, tri));
}

TEST(ElementTest, SubdivideParallelogram)
{
    const Parallelogram quad(Node(-2, -1), Node(1, -1), Node(2, 1), Node(-1, 1));

    const auto subdivision1 = quad.subdivide(Vector2mpq{0, 0});
    ASSERT_EQ(subdivision1.size(), 4);
    EXPECT_EQ(*subdivision1.at(0), Parallelogram(Node("-1/2", "-1"), Node(0, 0), Node("-3/2", "0"), Node(-2, -1)));
    EXPECT_EQ(*subdivision1.at(1), Parallelogram(Node("3/2", "0"), Node(0, 0), Node("-1/2", "-1"), Node(1, -1)));
    EXPECT_EQ(*subdivision1.at(2), Parallelogram(Node("1/2", "1"), Node(0, 0), Node("3/2", "0"), Node(2, 1)));
    EXPECT_EQ(*subdivision1.at(3), Parallelogram(Node("-3/2", "0"), Node(0, 0), Node("1/2", "1"), Node(-1, 1)));

    const auto subdivision2 = quad.subdivide(Vector2mpq{"1/2", "-1"});
    ASSERT_EQ(subdivision2.size(), 2);
    EXPECT_EQ(*subdivision2.at(0), Parallelogram(Node(-2, -1), Node("1/2", "-1"), Node("3/2", "1"), Node(-1, 1)));
    EXPECT_EQ(*subdivision2.at(1), Parallelogram(Node("3/2", "1"), Node("1/2", "-1"), Node(1, -1), Node(2, 1)));

    const auto subdivision3 = quad.subdivide(Vector2mpq{"7/4", "1/2"});
    ASSERT_EQ(subdivision3.size(), 2);
    EXPECT_EQ(*subdivision3.at(0), Parallelogram(Node(1, -1), Node("7/4", "1/2"), Node("-5/4", "1/2"), Node(-2, -1)));
    EXPECT_EQ(*subdivision3.at(1), Parallelogram(Node("-5/4", "1/2"), Node("7/4", "1/2"), Node(2, 1), Node(-1, 1)));

    const auto subdivision4 = quad.subdivide(Vector2mpq{0, 1});
    ASSERT_EQ(subdivision4.size(), 2);
    EXPECT_EQ(*subdivision4.at(0), Parallelogram(Node(2, 1), Node(0, 1), Node(-1, -1), Node(1, -1)));
    EXPECT_EQ(*subdivision4.at(1), Parallelogram(Node(-1, -1), Node(0, 1), Node(-1, 1), Node(-2, -1)));

    const auto subdivision5 = quad.subdivide(Vector2mpq{"-3/2", "0"});
    ASSERT_EQ(subdivision5.size(), 2);
    EXPECT_EQ(*subdivision5.at(0), Parallelogram(Node(-1, 1), Node("-3/2", "0"), Node("3/2", "0"), Node(2, 1)));
    EXPECT_EQ(*subdivision5.at(1), Parallelogram(Node("3/2", "0"), Node("-3/2", "0"), Node(-2, -1), Node(1, -1)));

    const auto subdivision6 = quad.subdivide(Vector2mpq{2, 1});
    ASSERT_EQ(subdivision6.size(), 1);
    EXPECT_EQ(*subdivision6.at(0), quad);

    const auto subdivision7 = quad.subdivide(Vector2mpq{-1, 1});
    ASSERT_EQ(subdivision7.size(), 1);
    EXPECT_EQ(*subdivision7.at(0), quad);

    const auto subdivision8 = quad.subdivide(Vector2mpq{-2, -1});
    ASSERT_EQ(subdivision8.size(), 1);
    EXPECT_EQ(*subdivision8.at(0), quad);

    const auto subdivision9 = quad.subdivide(Vector2mpq{1, -1});
    ASSERT_EQ(subdivision9.size(), 1);
    EXPECT_EQ(*subdivision9.at(0), quad);
}

TEST(ElementTest, SubdivideTriangle)
{
    const Triangle tri(Node(0, 0), Node(1, 0), Node(0, 1));

    const auto subdivision1 = tri.subdivide(Vector2mpq{"1/2", "1/4"});
    ASSERT_EQ(subdivision1.size(), 3);
    EXPECT_EQ(*subdivision1.at(0), Triangle(Node(0, 0), Node("1/2", "1/4"), Node(0, 1)));
    EXPECT_EQ(*subdivision1.at(1), Triangle(Node(1, 0), Node("1/2", "1/4"), Node(0, 0)));
    EXPECT_EQ(*subdivision1.at(2), Triangle(Node(0, 1), Node("1/2", "1/4"), Node(1, 0)));

    const auto subdivision2 = tri.subdivide(Vector2mpq{"1/2", "1/2"});
    ASSERT_EQ(subdivision2.size(), 2);
    EXPECT_EQ(*subdivision2.at(0), Triangle(Node(1, 0), Node("1/2", "1/2"), Node(0, 0)));
    EXPECT_EQ(*subdivision2.at(1), Triangle(Node(0, 0), Node("1/2", "1/2"), Node(0, 1)));

    const auto subdivision3 = tri.subdivide(Vector2mpq{"0", "3/4"});
    ASSERT_EQ(subdivision3.size(), 2);
    EXPECT_EQ(*subdivision3.at(0), Triangle(Node(0, 1), Node("0", "3/4"), Node(1, 0)));
    EXPECT_EQ(*subdivision3.at(1), Triangle(Node(1, 0), Node("0", "3/4"), Node(0, 0)));

    const auto subdivision4 = tri.subdivide(Vector2mpq{"1/2", "0"});
    ASSERT_EQ(subdivision4.size(), 2);
    EXPECT_EQ(*subdivision4.at(0), Triangle(Node(0, 0), Node("1/2", "0"), Node(0, 1)));
    EXPECT_EQ(*subdivision4.at(1), Triangle(Node(0, 1), Node("1/2", "0"), Node(1, 0)));

    const auto subdivision5 = tri.subdivide(Vector2mpq{0, 0});
    ASSERT_EQ(subdivision5.size(), 1);
    EXPECT_EQ(*subdivision5.at(0), tri);

    const auto subdivision6 = tri.subdivide(Vector2mpq{1, 0});
    ASSERT_EQ(subdivision6.size(), 1);
    EXPECT_EQ(*subdivision6.at(0), tri);

    const auto subdivision7 = tri.subdivide(Vector2mpq{0, 1});
    ASSERT_EQ(subdivision7.size(), 1);
    EXPECT_EQ(*subdivision7.at(0), tri);
}

TEST(ElementTest, SubdivideAbortsIfPointOutsideElement)
{
    const Parallelogram quad(Node(-2, -1), Node(1, -1), Node(2, 1), Node(-1, 1));
    EXPECT_DEATH(quad.subdivide(Vector2mpq{3, 0}), ".*");

    const Triangle tri(Node(0, 0), Node(1, 0), Node(0, 1));
    EXPECT_DEATH(tri.subdivide(Vector2mpq{-1, 0}), ".*");
}
} // namespace fem::ut
