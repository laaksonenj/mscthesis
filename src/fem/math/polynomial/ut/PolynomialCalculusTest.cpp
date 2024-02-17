#include <gtest/gtest.h>

#include "fem/math/polynomial/PolynomialCalculus.hpp"

namespace fem::ut
{
TEST(CalculusTest, Differentiation)
{
    EXPECT_EQ(diff(Polynomial1D(4)), Polynomial1D(0));
    EXPECT_EQ(diff(Polynomial1D("-1/4t^2+1/4t-1")), Polynomial1D("-1/2t+1/4"));
    EXPECT_EQ(diff(Polynomial1D("3/2t^3-3/2t")), Polynomial1D("9/2t^2-3/2"));
    EXPECT_EQ(diff(Polynomial2D("9/2xy^2+3/2y^2-x-1"), 'x'), Polynomial2D("9/2y^2-1"));
    EXPECT_EQ(diff(Polynomial2D("9/2xy^2+3/2y^2-x-1"), 'y'), Polynomial2D("9xy+3y"));
}

TEST(CalculusTest, IntegratePolynomialOverReferenceElements)
{
    EXPECT_EQ(integrateOverReferenceElement(Polynomial2D("-1/4x^2+1/4x-1"), ElementType_Parallelogram), mpq_class("-13/3"));
    EXPECT_EQ(integrateOverReferenceElement(mpq_class("-3/8") * Polynomial2D("y^4+4xy^3+3y^3+6x^2y^2+9xy^2+2y^2-2x^2-2xy-3x-y-1"), ElementType_Parallelogram), mpq_class("1/5"));
    EXPECT_EQ(integrateOverReferenceElement(Polynomial2D(0), ElementType_Parallelogram), mpq_class(0));
    EXPECT_EQ(integrateOverReferenceElement(Polynomial2D(1), ElementType_Parallelogram), mpq_class(4));

    EXPECT_EQ(integrateOverReferenceElement(Polynomial2D("3/2y^3-3/2y"), ElementType_Triangle), mpq_class("-7/40"));
    EXPECT_EQ(integrateOverReferenceElement(-12 * Polynomial2D("6y^2+x^2+6xy-2x-6y+1"), ElementType_Triangle), mpq_class(0));
    EXPECT_EQ(integrateOverReferenceElement(Polynomial2D(0), ElementType_Triangle), mpq_class(0));
    EXPECT_EQ(integrateOverReferenceElement(Polynomial2D(1), ElementType_Triangle), mpq_class("1/2"));
}
} // namespace fem::ut
