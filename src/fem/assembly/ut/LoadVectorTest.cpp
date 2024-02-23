#include <gtest/gtest.h>

#include "fem/assembly/LoadVector.hpp"

namespace fem::ut
{
TEST(LoadVectorTest, ExtractSubLoadVector)
{
    const std::vector<Node> nodes = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, 2}};
    const std::vector<std::vector<Mesh::ElementIndex>> elements = {{0, 1, 2, 3}, {3, 2, 4}};
    const Mesh mesh(nodes, elements);
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Product;
    const Vector2mpq x_0{0, 0};
    for (int p_max = 1; p_max <= 4; p_max++)
    {
        const BasisFunctionIndexer superBasisFunctionIndexer(mesh, p_max, polynomialSpaceType);
        const VectorXmpq loadVector = assembleDiracLoadVector(superBasisFunctionIndexer, x_0);
        for (int p = 1; p <= p_max; p++)
        {
            const BasisFunctionIndexer subBasisFunctionIndexer(mesh, p, polynomialSpaceType);
            const VectorXmpq subLoadVector = assembleDiracLoadVector(subBasisFunctionIndexer, x_0);
            EXPECT_EQ(extractSubLoadVector(loadVector, superBasisFunctionIndexer, subBasisFunctionIndexer), subLoadVector);
        }
    }
}
} // namespace fem::ut
