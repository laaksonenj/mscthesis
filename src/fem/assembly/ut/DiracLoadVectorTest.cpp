#include <gtest/gtest.h>

#include "fem/assembly/DiracLoadVector.hpp"
#include "fem/assembly/ut/refdata/dirac_load1/RefDiracLoadVector1.hpp"
#include "fem/assembly/ut/refdata/dirac_load2/RefDiracLoadVector2.hpp"
#include "fem/assembly/ut/refdata/dirac_load3/RefDiracLoadVector3.hpp"

namespace fem::ut
{
TEST(DiracLoadVectorTest, DiracLoadVector1)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/dirac_load1/mesh.txt"};
    const uint32_t p = 4;
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Trunk;
    const FemContext ctx(std::make_shared<Mesh>(createMeshFromFile(meshFilename)), p, polynomialSpaceType);
    const Vector2mpq x_0{1, 2};
    ShapeFunctionFactory shapeFunctionFactory;
    shapeFunctionFactory.createShapeFunctions(ElementType_Parallelogram, p);
    shapeFunctionFactory.createShapeFunctions(ElementType_Triangle, p);
    EXPECT_EQ(assembleDiracLoadVector(ctx, x_0), refdata::refDiracLoadVector1);
    EXPECT_EQ(assembleDiracLoadVector(ctx, x_0, shapeFunctionFactory), refdata::refDiracLoadVector1);
}

TEST(DiracLoadVectorTest, DiracLoadVector2)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/dirac_load2/mesh.txt"};
    const uint32_t p = 4;
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Product;
    const FemContext ctx(std::make_shared<Mesh>(createMeshFromFile(meshFilename)), p, polynomialSpaceType);
    const Vector2mpq x_0{mpq_class("3/2"), mpq_class("5/2")};
    ShapeFunctionFactory shapeFunctionFactory;
    shapeFunctionFactory.createShapeFunctions(ElementType_Parallelogram, p);
    shapeFunctionFactory.createShapeFunctions(ElementType_Triangle, p);
    EXPECT_EQ(assembleDiracLoadVector(ctx, x_0), refdata::refDiracLoadVector2);
    EXPECT_EQ(assembleDiracLoadVector(ctx, x_0, shapeFunctionFactory), refdata::refDiracLoadVector2);
}

TEST(DiracLoadVectorTest, DiracLoadVector3)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/dirac_load3/mesh.txt"};
    const uint32_t p = 3;
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Trunk;
    const FemContext ctx(std::make_shared<Mesh>(createMeshFromFile(meshFilename)), p, polynomialSpaceType);
    const Vector2mpq x_0{mpq_class(1), mpq_class("5/2")};
    ShapeFunctionFactory shapeFunctionFactory;
    shapeFunctionFactory.createShapeFunctions(ElementType_Parallelogram, p);
    shapeFunctionFactory.createShapeFunctions(ElementType_Triangle, p);
    EXPECT_EQ(assembleDiracLoadVector(ctx, x_0), refdata::refDiracLoadVector3);
    EXPECT_EQ(assembleDiracLoadVector(ctx, x_0, shapeFunctionFactory), refdata::refDiracLoadVector3);
}
} // namespace fem::ut
