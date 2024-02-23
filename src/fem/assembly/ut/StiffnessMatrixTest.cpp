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
    const BasisFunctionIndexer basisFunctionIndexer(mesh, p, polynomialSpaceType);
    const ShapeFunctionFactory shapeFunctionFactory;
    EXPECT_EQ(assembleStiffnessMatrix(basisFunctionIndexer), refdata::refStiffnessMatrix1);
    EXPECT_EQ(assembleStiffnessMatrix(basisFunctionIndexer, shapeFunctionFactory), refdata::refStiffnessMatrix1);
}

TEST(StiffnessMatrixTest, StiffnessMatrix2)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/stiffness_matrix2/mesh.txt"};
    const Mesh mesh = createMeshFromFile(meshFilename);
    const uint32_t p = 4;
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Product;
    const BasisFunctionIndexer basisFunctionIndexer(mesh, p, polynomialSpaceType);
    const ShapeFunctionFactory shapeFunctionFactory;
    EXPECT_EQ(assembleStiffnessMatrix(basisFunctionIndexer), refdata::refStiffnessMatrix2);
    EXPECT_EQ(assembleStiffnessMatrix(basisFunctionIndexer, shapeFunctionFactory), refdata::refStiffnessMatrix2);
}

TEST(StiffnessMatrixTest, StiffnessMatrix3)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/stiffness_matrix3/mesh.txt"};
    const Mesh mesh = createMeshFromFile(meshFilename);
    const uint32_t p = 1;
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Trunk;
    const BasisFunctionIndexer basisFunctionIndexer(mesh, p, polynomialSpaceType);
    const ShapeFunctionFactory shapeFunctionFactory;
    EXPECT_EQ(assembleStiffnessMatrix(basisFunctionIndexer), refdata::refStiffnessMatrix3);
    EXPECT_EQ(assembleStiffnessMatrix(basisFunctionIndexer, shapeFunctionFactory), refdata::refStiffnessMatrix3);
}

TEST(StiffnessMatrixTest, ExtractSubStiffnessMatrix)
{
    const std::vector<Node> nodes = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, 2}};
    const std::vector<std::vector<Mesh::ElementIndex>> elements = {{0, 1, 2, 3}, {3, 2, 4}};
    const Mesh mesh(nodes, elements);
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Product;
    for (int p_max = 1; p_max <= 4; p_max++)
    {
        const BasisFunctionIndexer superBasisFunctionIndexer(mesh, p_max, polynomialSpaceType);
        const MatrixXmpq stiffnessMatrix = assembleStiffnessMatrix(superBasisFunctionIndexer);
        for (int p = 1; p <= p_max; p++)
        {
            const BasisFunctionIndexer subBasisFunctionIndexer(mesh, p, polynomialSpaceType);
            const MatrixXmpq subStiffnessMatrix = assembleStiffnessMatrix(subBasisFunctionIndexer);
            EXPECT_EQ(extractSubStiffnessMatrix(stiffnessMatrix, superBasisFunctionIndexer, subBasisFunctionIndexer), subStiffnessMatrix);
        }
    }
}
} // namespace fem::ut
