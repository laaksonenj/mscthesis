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
    const FemContext ctx(std::make_shared<Mesh>(createMeshFromFile(meshFilename)), 4, PolynomialSpaceType_Trunk);
    ShapeFunctionFactory shapeFunctionFactory;
    EXPECT_EQ(assembleStiffnessMatrix(ctx), refdata::refStiffnessMatrix1);
    EXPECT_EQ(assembleStiffnessMatrix(ctx, shapeFunctionFactory), refdata::refStiffnessMatrix1);
}

TEST(StiffnessMatrixTest, StiffnessMatrix2)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/stiffness_matrix2/mesh.txt"};
    const FemContext ctx(std::make_shared<Mesh>(createMeshFromFile(meshFilename)), 4, PolynomialSpaceType_Product);
    ShapeFunctionFactory shapeFunctionFactory;
    EXPECT_EQ(assembleStiffnessMatrix(ctx), refdata::refStiffnessMatrix2);
    EXPECT_EQ(assembleStiffnessMatrix(ctx, shapeFunctionFactory), refdata::refStiffnessMatrix2);
}

TEST(StiffnessMatrixTest, StiffnessMatrix3)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/stiffness_matrix3/mesh.txt"};
    const FemContext ctx(std::make_shared<Mesh>(createMeshFromFile(meshFilename)), 1, PolynomialSpaceType_Trunk);
    ShapeFunctionFactory shapeFunctionFactory;
    EXPECT_EQ(assembleStiffnessMatrix(ctx), refdata::refStiffnessMatrix3);
    EXPECT_EQ(assembleStiffnessMatrix(ctx, shapeFunctionFactory), refdata::refStiffnessMatrix3);
}

TEST(StiffnessMatrixTest, ExtractSubStiffnessMatrix)
{
    const std::vector<Node> nodes = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, 2}};
    const std::vector<std::vector<Mesh::ElementIndex>> elements = {{0, 1, 2, 3}, {3, 2, 4}};
    const auto mesh = std::make_shared<Mesh>(nodes, elements);
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Product;
    for (int p_max = 1; p_max <= 4; p_max++)
    {
        const FemContext ctx(mesh, p_max, polynomialSpaceType);
        const MatrixXmpq stiffnessMatrix = assembleStiffnessMatrix(ctx);
        for (int p = 1; p <= p_max; p++)
        {
            const FemContext subCtx(mesh, p, polynomialSpaceType);
            const MatrixXmpq subStiffnessMatrix = assembleStiffnessMatrix(subCtx);
            EXPECT_EQ(extractSubStiffnessMatrix(stiffnessMatrix, ctx, p), subStiffnessMatrix);
        }
    }
}
} // namespace fem::ut
