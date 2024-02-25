#include <gtest/gtest.h>

#include "fem/assembly/NeumannLoadVector.hpp"
#include "fem/assembly/ut/ReferenceBoundaryFunctions.hpp"
#include "fem/assembly/ut/refdata/neumann_load1/RefNeumannLoadVector1.hpp"
#include "fem/assembly/ut/refdata/neumann_load2/RefNeumannLoadVector2.hpp"

namespace fem::ut
{
namespace
{
void checkNearElementwise(const VectorXmpq& lhs, const VectorXmpq& rhs)
{
    EXPECT_EQ(lhs.size(), rhs.size());
    for (int i = 0; i < lhs.size(); i++)
    {
        EXPECT_NEAR(lhs(i).get_d(), rhs(i).get_d(), 1e-7);
    }
}
} // namespace

TEST(NeumannLoadVectorTest, NeumannLoadVector1)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/neumann_load1/mesh.txt"};
    const Mesh mesh = createMeshFromFile(meshFilename);
    const uint32_t p = 4;
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Trunk;
    const BivariateFunction& g = referenceBoundaryFunctions.at(0);
    const uint32_t elementIdx = 1;
    const uint32_t localSideIdx = 0;
    const BasisFunctionIndexer basisFunctionIndexer(mesh, p, polynomialSpaceType);
    const ShapeFunctionFactory shapeFunctionFactory;
    checkNearElementwise(assembleNeumannLoadVector(basisFunctionIndexer, g, elementIdx, localSideIdx), refdata::refNeumannLoadVector1);
    checkNearElementwise(assembleNeumannLoadVector(basisFunctionIndexer, g, elementIdx, localSideIdx, shapeFunctionFactory), refdata::refNeumannLoadVector1);
}

TEST(NeumannLoadVectorTest, NeumannLoadVector2)
{
    const std::string meshFilename = std::string{SRC_DIR} + std::string{"/refdata/neumann_load2/mesh.txt"};
    const Mesh mesh = createMeshFromFile(meshFilename);
    const uint32_t p = 4;
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Product;
    const BivariateFunction& g = referenceBoundaryFunctions.at(1);
    const uint32_t elementIdx = 1;
    const uint32_t localSideIdx = 3;
    const BasisFunctionIndexer basisFunctionIndexer(mesh, p, polynomialSpaceType);
    const ShapeFunctionFactory shapeFunctionFactory;
    checkNearElementwise(assembleNeumannLoadVector(basisFunctionIndexer, g, elementIdx, localSideIdx), refdata::refNeumannLoadVector2);
    checkNearElementwise(assembleNeumannLoadVector(basisFunctionIndexer, g, elementIdx, localSideIdx, shapeFunctionFactory), refdata::refNeumannLoadVector2);
}
} // namespace fem::ut