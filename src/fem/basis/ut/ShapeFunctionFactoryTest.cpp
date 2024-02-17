#include <gtest/gtest.h>

#include "fem/basis/ShapeFunctionFactory.hpp"

namespace fem::ut
{
namespace
{
ShapeFunctionFactory f;
} // namespace

TEST(ShapeFunctionFactoryTest, NodalShapeFunctionsTriangle)
{
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, NodalShapeFunctionDescriptor(0)), Polynomial2D("1-x-y"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, NodalShapeFunctionDescriptor(1)), Polynomial2D("x"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, NodalShapeFunctionDescriptor(2)), Polynomial2D("y"));
}

TEST(ShapeFunctionFactoryTest, SideShapeFunctionsTriangle)
{
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, SideShapeFunctionDescriptor(0, 2)), Polynomial2D("4 x^2 + 4 x y - 4 x"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, SideShapeFunctionDescriptor(0, 3)), Polynomial2D("24 x^3 + 36 x^2 y - 36 x^2 + 12 x y^2 - 24 x y + 12 x"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, SideShapeFunctionDescriptor(0, 4)), Polynomial2D("120 x^4 + 240 x^3 y - 240 x^3 + 150 x^2 y^2 - 300 x^2 y + 144 x^2 + 30 x y^3 - 90 x y^2 + 84 x y - 24 x"));

    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, SideShapeFunctionDescriptor(1, 2)), Polynomial2D("-4xy"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, SideShapeFunctionDescriptor(1, 3)), Polynomial2D("12 x^2 y - 12 x y^2"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, SideShapeFunctionDescriptor(1, 4)), Polynomial2D("-30 x^3 y + 60 x^2 y^2 - 30 x y^3 + 6 x y"));

    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, SideShapeFunctionDescriptor(2, 2)), Polynomial2D("4 x y + 4 y^2 - 4 y"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, SideShapeFunctionDescriptor(2, 3)), Polynomial2D("-12 x^2 y - 36 x y^2 + 24 x y - 24 y^3 + 36 y^2 - 12 y"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, SideShapeFunctionDescriptor(2, 4)), Polynomial2D("30 x^3 y + 150 x^2 y^2 - 90 x^2 y + 240 x y^3 - 300 x y^2 + 84 x y + 120 y^4 - 240 y^3 + 144 y^2 - 24 y"));
}

TEST(ShapeFunctionFactoryTest, InternalShapeFunctionsTriangle)
{
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, InternalShapeFunctionDescriptor(0, 0)), Polynomial2D("-x^2 y - x y^2 + x y"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, InternalShapeFunctionDescriptor(0, 1)), Polynomial2D("-2 x^2 y^2 + x^2 y - 2 x y^3 + 3 x y^2 - x y"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, InternalShapeFunctionDescriptor(1, 0)), Polynomial2D("-2 x^3 y - 2 x^2 y^2 + 3 x^2 y + x y^2 - x y"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, InternalShapeFunctionDescriptor(1, 1)), Polynomial2D("-4 x^3 y^2 + 2 x^3 y - 4 x^2 y^3 + 8 x^2 y^2 - 3 x^2 y + 2 x y^3 - 3 x y^2 + x y"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Triangle, InternalShapeFunctionDescriptor(2, 3)), Polynomial2D("-120 x^4 y^4 + 180 x^4 y^3 - 72 x^4 y^2 + 6 x^4 y - 120 x^3 y^5 + 420 x^3 y^4 - 432 x^3 y^3 + 150 x^3 y^2 - 12 x^3 y + 120 x^2 y^5 - 320 x^2 y^4 + 282 x^2 y^3 - 90 x^2 y^2 + 7 x^2 y - 20 x y^5 + 50 x y^4 - 42 x y^3 + 13 x y^2 - x y"));
}

TEST(ShapeFunctionFactoryTest, NodalShapeFunctionsQuadrilateral)
{
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, NodalShapeFunctionDescriptor(0)), Polynomial2D("1/4xy-1/4x-1/4y+1/4"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, NodalShapeFunctionDescriptor(1)), Polynomial2D("-1/4xy+1/4x-1/4y+1/4"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, NodalShapeFunctionDescriptor(2)), Polynomial2D("1/4xy+1/4x+1/4y+1/4"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, NodalShapeFunctionDescriptor(3)), Polynomial2D("-1/4xy-1/4x+1/4y+1/4"));
}

TEST(ShapeFunctionFactoryTest, SideShapeFunctionsQuadrilateral)
{
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(0, 2)), Polynomial2D("-1/2x^2y+1/2x^2+1/2y-1/2"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(0, 3)), Polynomial2D("-3/2x^3y+3/2x^3+3/2xy-3/2x"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(0, 4)), Polynomial2D("-15/4x^4y+15/4x^4+9/2x^2y-9/2x^2-3/4y+3/4"));

    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(1, 2)), Polynomial2D("1/2xy^2-1/2x+1/2y^2-1/2"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(1, 3)), Polynomial2D("3/2xy^3-3/2xy+3/2y^3-3/2y"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(1, 4)), Polynomial2D("15/4xy^4-9/2xy^2+3/4x+15/4y^4-9/2y^2+3/4"));

    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(2, 2)), Polynomial2D("1/2x^2y+1/2x^2-1/2y-1/2"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(2, 3)), Polynomial2D("-3/2x^3y-3/2x^3+3/2xy+3/2x"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(2, 4)), Polynomial2D("15/4x^4y+15/4x^4-9/2x^2y-9/2x^2+3/4y+3/4"));

    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(3, 2)), Polynomial2D("-1/2xy^2+1/2x+1/2y^2-1/2"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(3, 3)), Polynomial2D("3/2xy^3-3/2xy-3/2y^3+3/2y"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(3, 4)), Polynomial2D("-15/4xy^4+9/2xy^2-3/4x+15/4y^4-9/2y^2+3/4"));
}

TEST(ShapeFunctionFactoryTest, InternalShapeFunctionsQuadrilateral)
{
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, InternalShapeFunctionDescriptor(2, 2)), Polynomial2D("x^2 y^2 - x^2 - y^2 + 1"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, InternalShapeFunctionDescriptor(2, 3)), Polynomial2D("3 x^2 y^3 - 3 x^2 y - 3 y^3 + 3 y"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, InternalShapeFunctionDescriptor(3, 2)), Polynomial2D("3 x^3 y^2 - 3 x^3 - 3 x y^2 + 3 x"));
    EXPECT_EQ(f.getShapeFunction(ElementType_Parallelogram, InternalShapeFunctionDescriptor(3, 3)), Polynomial2D("9 x^3 y^3 - 9 x^3 y - 9 x y^3 + 9 x y"));
}

TEST(ShapeFunctionFactoryTest, ShapeFunctionDerivatives)
{
    const Polynomial2D& p1 = f.getShapeFunction(ElementType_Triangle, NodalShapeFunctionDescriptor(0));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Triangle, NodalShapeFunctionDescriptor(0), 'x'), diff(p1, 'x'));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Triangle, NodalShapeFunctionDescriptor(0), 'y'), diff(p1, 'y'));

    const Polynomial2D& p2 = f.getShapeFunction(ElementType_Triangle, SideShapeFunctionDescriptor(1, 2));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Triangle, SideShapeFunctionDescriptor(1, 2), 'x'), diff(p2, 'x'));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Triangle, SideShapeFunctionDescriptor(1, 2), 'y'), diff(p2, 'y'));

    const Polynomial2D& p3 = f.getShapeFunction(ElementType_Triangle, InternalShapeFunctionDescriptor(2, 3));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Triangle, InternalShapeFunctionDescriptor(2, 3), 'x'), diff(p3, 'x'));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Triangle, InternalShapeFunctionDescriptor(2, 3), 'y'), diff(p3, 'y'));

    const Polynomial2D& p4 = f.getShapeFunction(ElementType_Parallelogram, NodalShapeFunctionDescriptor(2));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Parallelogram, NodalShapeFunctionDescriptor(2), 'x'), diff(p4, 'x'));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Parallelogram, NodalShapeFunctionDescriptor(2), 'y'), diff(p4, 'y'));

    const Polynomial2D& p5 = f.getShapeFunction(ElementType_Parallelogram, SideShapeFunctionDescriptor(0, 4));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Parallelogram, SideShapeFunctionDescriptor(0, 4), 'x'), diff(p5, 'x'));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Parallelogram, SideShapeFunctionDescriptor(0, 4), 'y'), diff(p5, 'y'));

    const Polynomial2D& p6 = f.getShapeFunction(ElementType_Parallelogram, InternalShapeFunctionDescriptor(2, 2));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Parallelogram, InternalShapeFunctionDescriptor(2, 2), 'x'), diff(p6, 'x'));
    EXPECT_EQ(f.getShapeFunctionDerivative(ElementType_Parallelogram, InternalShapeFunctionDescriptor(2, 2), 'y'), diff(p6, 'y'));
}
} // fem::ut
