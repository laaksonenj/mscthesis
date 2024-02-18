#include <gtest/gtest.h>

#include "fem/assembly/StiffnessMatrix.hpp"
#include "fem/assembly/ut/refdata/stiffness_matrix1/RefStiffnessMatrix1.hpp"
#include "fem/assembly/ut/refdata/stiffness_matrix2/RefStiffnessMatrix2.hpp"
#include "fem/assembly/ut/refdata/stiffness_matrix3/RefStiffnessMatrix3.hpp"

namespace fem::ut
{
TEST(StiffnessMatrixTest, StiffnessMatrix1)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/stiffness_matrix1/mesh.txt"};
    const Mesh mesh = createMeshFromFile(meshFilename);
    const uint32_t p = 4;
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Trunk;
    const ShapeFunctionFactory shapeFunctionFactory;
    EXPECT_EQ(assembleStiffnessMatrix(mesh, p, polynomialSpaceType), refdata::refStiffnessMatrix1);
    EXPECT_EQ(assembleStiffnessMatrix(mesh, p, polynomialSpaceType, shapeFunctionFactory), refdata::refStiffnessMatrix1);
}

TEST(StiffnessMatrixTest, StiffnessMatrix2)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/stiffness_matrix2/mesh.txt"};
    const Mesh mesh = createMeshFromFile(meshFilename);
    const uint32_t p = 4;
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Product;
    const ShapeFunctionFactory shapeFunctionFactory;
    EXPECT_EQ(assembleStiffnessMatrix(mesh, p, polynomialSpaceType), refdata::refStiffnessMatrix2);
    EXPECT_EQ(assembleStiffnessMatrix(mesh, p, polynomialSpaceType, shapeFunctionFactory), refdata::refStiffnessMatrix2);
}

TEST(StiffnessMatrixTest, StiffnessMatrix3)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/stiffness_matrix3/mesh.txt"};
    const Mesh mesh = createMeshFromFile(meshFilename);
    const uint32_t p = 1;
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Trunk;
    const ShapeFunctionFactory shapeFunctionFactory;
    EXPECT_EQ(assembleStiffnessMatrix(mesh, p, polynomialSpaceType), refdata::refStiffnessMatrix3);
    EXPECT_EQ(assembleStiffnessMatrix(mesh, p, polynomialSpaceType, shapeFunctionFactory), refdata::refStiffnessMatrix3);
}
} // namespace fem::ut
