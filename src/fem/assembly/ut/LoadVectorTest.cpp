#include <gtest/gtest.h>

#include "fem/assembly/LoadVector.hpp"

namespace fem::ut
{
TEST(LoadVectorTest, ExtractSubLoadVector)
{
    const std::vector<Node> nodes = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, 2}};
    const std::vector<std::vector<Mesh::ElementIndex>> elements = {{0, 1, 2, 3}, {3, 2, 4}};
    const auto mesh = std::make_shared<Mesh>(nodes, elements);
    const PolynomialSpaceType polynomialSpaceType = PolynomialSpaceType_Product;
    const Vector2mpq x_0{0, 0};
    for (int p_max = 1; p_max <= 4; p_max++)
    {
        const FemContext ctx(mesh, p_max, polynomialSpaceType);
        const VectorXmpq loadVector = assembleDiracLoadVector(ctx, x_0);
        for (int p = 1; p <= p_max; p++)
        {
            const FemContext subCtx(mesh, p, polynomialSpaceType);
            const VectorXmpq subLoadVector = assembleDiracLoadVector(subCtx, x_0);
            EXPECT_EQ(extractSubLoadVector(ctx, loadVector, p), subLoadVector);
        }
    }
}
} // namespace fem::ut
