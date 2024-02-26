#include <gtest/gtest.h>

#include "fem/basis/BasisFunctionFactory.hpp"

namespace fem::ut
{
class BasisFunctionFactoryTest : public testing::Test
{
protected:
    const std::vector<Node> vertices{
        {"0", "-1"},
        {"1", "-1"},
        {"1", "0"},
        {"0", "1"},
        {"-1", "1"},
        {"-1", "0"},
        {"0", "0"}
    };
    const std::vector<std::vector<Mesh::NodeIndex>> elements{
        {0, 1, 2, 6},
        {2, 3, 6},
        {3, 4, 5, 6},
        {6, 5, 0}
    };
    std::shared_ptr<Mesh> mesh{new Mesh(vertices, elements)};

    void checkContinuityOverSides(uint32_t p)
    {
        const FemContext ctx(mesh, p, PolynomialSpaceType_Product);
        BasisFunctionFactory factory(ctx);
        checkContinuityOverSide2(factory, p);
        checkContinuityOverSide3(factory, p);
        checkContinuityOverSide5(factory, p);
        checkContinuityOverSide8(factory, p);
    }

    void checkContinuityOverSide2(const BasisFunctionFactory& factory, uint32_t p)
    {
        const uint32_t elementIdx1 = 0;
        const uint32_t elementIdx2 = 1;
        const Polynomial1D sideParameterizationX("t");
        const Polynomial1D sideParameterizationY(0);
        const AffineMap G1 = mesh->getElement(elementIdx1).getReferenceElementMap().inverse();
        const AffineMap G2 = mesh->getElement(elementIdx2).getReferenceElementMap().inverse();
        
        /* Nodal basis functions */
        EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 2), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(factory.getShapeFunction(elementIdx2, 0), G2), sideParameterizationX, sideParameterizationY));
        EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 3), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(factory.getShapeFunction(elementIdx2, 2), G2), sideParameterizationX, sideParameterizationY));

        /* Side basis functions */
        for (int i = 0; i < p-1; i++)
        {
            EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 4 + 2*(p-1) + i), G1), sideParameterizationX, sideParameterizationY),
                      compose(compose(factory.getShapeFunction(elementIdx2, 3 + 2*(p-1) + i), G2), sideParameterizationX, sideParameterizationY));
        }
    }
    
    void checkContinuityOverSide3(const BasisFunctionFactory& factory, uint32_t p)
    {
        const uint32_t elementIdx1 = 0;
        const uint32_t elementIdx2 = 3;
        const Polynomial1D sideParameterizationX(0);
        const Polynomial1D sideParameterizationY("t-1");
        const AffineMap G1 = mesh->getElement(elementIdx1).getReferenceElementMap().inverse();
        const AffineMap G2 = mesh->getElement(elementIdx2).getReferenceElementMap().inverse();
        
        /* Nodal basis functions */
        EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 0), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(factory.getShapeFunction(elementIdx2, 2), G2), sideParameterizationX, sideParameterizationY));
        EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 3), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(factory.getShapeFunction(elementIdx2, 0), G2), sideParameterizationX, sideParameterizationY));

        /* Side basis functions */
        for (int i = 0; i < p-1; i++)
        {
            EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 4 + 3*(p-1) + i), G1), sideParameterizationX, sideParameterizationY),
                      compose(compose(factory.getShapeFunction(elementIdx2, 3 + 2*(p-1) + i), G2), sideParameterizationX, sideParameterizationY));
        }
    }

    void checkContinuityOverSide5(const BasisFunctionFactory& factory, uint32_t p)
    {
        const uint32_t elementIdx1 = 1;
        const uint32_t elementIdx2 = 2;
        const Polynomial1D sideParameterizationX(0);
        const Polynomial1D sideParameterizationY("t");
        const AffineMap G1 = mesh->getElement(elementIdx1).getReferenceElementMap().inverse();
        const AffineMap G2 = mesh->getElement(elementIdx2).getReferenceElementMap().inverse();
        
        /* Nodal basis functions */
        EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 1), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(factory.getShapeFunction(elementIdx2, 0), G2), sideParameterizationX, sideParameterizationY));
        EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 2), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(factory.getShapeFunction(elementIdx2, 3), G2), sideParameterizationX, sideParameterizationY));

        /* Side basis functions */
        for (int i = 0; i < p-1; i++)
        {
            EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 3 + 1*(p-1) + i), G1), sideParameterizationX, sideParameterizationY),
                      compose(compose(factory.getShapeFunction(elementIdx2, 4 + 3*(p-1) + i), G2), sideParameterizationX, sideParameterizationY));
        }
    }

    void checkContinuityOverSide8(const BasisFunctionFactory& factory, uint32_t p)
    {
        const uint32_t elementIdx1 = 2;
        const uint32_t elementIdx2 = 3;
        const Polynomial1D sideParameterizationX("-t");
        const Polynomial1D sideParameterizationY(0);
        const AffineMap G1 = mesh->getElement(elementIdx1).getReferenceElementMap().inverse();
        const AffineMap G2 = mesh->getElement(elementIdx2).getReferenceElementMap().inverse();
        
        /* Nodal basis functions */
        EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 2), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(factory.getShapeFunction(elementIdx2, 1), G2), sideParameterizationX, sideParameterizationY));
        EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 3), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(factory.getShapeFunction(elementIdx2, 0), G2), sideParameterizationX, sideParameterizationY));

        /* Side basis functions */
        for (int i = 0; i < p-1; i++)
        {
            EXPECT_EQ(compose(compose(factory.getShapeFunction(elementIdx1, 4 + 2*(p-1) + i), G1), sideParameterizationX, sideParameterizationY),
                      compose(compose(factory.getShapeFunction(elementIdx2, 3 + 0*(p-1) + i), G2), sideParameterizationX, sideParameterizationY));
        }
    }
};

TEST_F(BasisFunctionFactoryTest, ContinuityOverSides)
{
    for (int p = 1; p <= 7; p++)
    {
        checkContinuityOverSides(p);
    }
}
} // namespace fem::ut
