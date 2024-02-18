#include <gtest/gtest.h>

#include "fem/assembly/StiffnessMatrix.hpp"

namespace fem::ut
{
TEST(StiffnessMatrixTest, StiffnessMatrixAssembly)
{
    const std::vector<Node> nodes{
        {"1", "0"},
        {"2", "1"},
        {"2", "3"},
        {"0", "3"},
        {"1", "2"}
    };
    const std::vector<std::vector<Mesh::NodeIndex>> elements{
        {2, 3, 4},
        {1, 2, 4, 0}
    };
    const Mesh mesh(nodes, elements);

    const uint32_t p = 4;
    const MatrixXmpq stiffnessMatrix = assembleStiffnessMatrix(mesh, p, PolynomialSpaceType_Trunk);
    EXPECT_EQ(stiffnessMatrix.rows(), 27);
    EXPECT_EQ(stiffnessMatrix.cols(), 27);

    /* Edge 2, k=3 vs. node 2 */
    uint32_t i = 5 + 2*(p-1) + 1;
    uint32_t j = 2;
    EXPECT_EQ(stiffnessMatrix(i, j), mpq_class("1/5"));
    EXPECT_EQ(stiffnessMatrix(j, i), mpq_class("1/5"));
}
} // namespace fem::ut
