#include <gtest/gtest.h>

#include <set>

#include "fem/basis/BasisFunctionIndexer.hpp"

namespace fem::ut
{
class BasisFunctionIndexerTest : public testing::Test
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
    const std::vector<std::vector<Mesh::SideIndex>> globalSideIndices{
        {0, 1, 2, 3},
        {4, 5, 2},
        {6, 7, 8, 5},
        {8, 9, 3}
    };
    const Mesh mesh{vertices, elements};

    std::set<uint32_t> basisFunctionIndices;

    void checkBasisFunctions(uint32_t p, PolynomialSpaceType polynomialSpaceType)
    {
        BasisFunctionIndexer indexer(mesh, p, polynomialSpaceType);
        EXPECT_EQ(indexer.getNumOfBasisFunctions(), getNumOfBasisFunctions(p, polynomialSpaceType));
        checkNodalBasisFunctions(indexer);
        checkSideBasisFunctions(indexer, p);
        checkInternalBasisFunctions(indexer, p, polynomialSpaceType);
        EXPECT_EQ(basisFunctionIndices.size(), indexer.getNumOfBasisFunctions());
        EXPECT_EQ(*basisFunctionIndices.rbegin(), indexer.getNumOfBasisFunctions() - 1);
    }

    uint32_t getNumOfBasisFunctions(uint32_t p, PolynomialSpaceType polynomialSpaceType)
    {
        uint32_t res = vertices.size();
        res += 10 * (p-1); // sides
        if (polynomialSpaceType == PolynomialSpaceType_Product)
        {
            res += 4 * (p-1)*(p-1);
        }
        else
        {
            res += 2 * (p-1)*(p-2)/2; // tris
            if (p >= 4)
            {
                res += 2 * (p-2)*(p-3)/2; // quads
            }
        }
        return res;
    }

    void checkNodalBasisFunctions(const BasisFunctionIndexer& indexer)
    {
        for (int elementIdx = 0; elementIdx < 4; elementIdx++)
        {
            for (int localVertexIdx = 0; localVertexIdx < elements.at(elementIdx).size(); localVertexIdx++)
            {
                BasisFunctionDescriptor desc = indexer.getBasisFunctionDescriptor(elementIdx, localVertexIdx);
                uint32_t basisFunctionIdx = indexer.getBasisFunctionIndex(elementIdx, localVertexIdx);
                EXPECT_EQ(desc, BasisFunctionDescriptor(NodalBasisFunctionDescriptor(elements.at(elementIdx).at(localVertexIdx))));
                EXPECT_EQ(basisFunctionIdx, indexer.getBasisFunctionIndex(desc));
                EXPECT_EQ(desc, indexer.getBasisFunctionDescriptor(basisFunctionIdx));
                basisFunctionIndices.insert(basisFunctionIdx);
            }
        }
    }

    void checkSideBasisFunctions(const BasisFunctionIndexer& indexer, uint32_t p)
    {
        for (int elementIdx = 0; elementIdx < 4; elementIdx++)
        {
            for (int localSideIdx = 0; localSideIdx < globalSideIndices.at(elementIdx).size(); localSideIdx++)
            {
                for (int k = 2; k <= p; k++)
                {
                    uint32_t shapeFunctionIdx = elements.at(elementIdx).size() + localSideIdx * (p-1) + k-2;
                    BasisFunctionDescriptor desc = indexer.getBasisFunctionDescriptor(elementIdx, shapeFunctionIdx);
                    uint32_t basisFunctionIdx = indexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
                    EXPECT_EQ(desc, BasisFunctionDescriptor(SideBasisFunctionDescriptor(globalSideIndices.at(elementIdx).at(localSideIdx), k)));
                    EXPECT_EQ(basisFunctionIdx, indexer.getBasisFunctionIndex(desc));
                    EXPECT_EQ(desc, indexer.getBasisFunctionDescriptor(basisFunctionIdx));
                    basisFunctionIndices.insert(basisFunctionIdx);
                }
            }
        }
    }

    void checkInternalBasisFunctions(const BasisFunctionIndexer& indexer, uint32_t p, PolynomialSpaceType polynomialSpaceType)
    {
        if (polynomialSpaceType == PolynomialSpaceType_Product)
        {
            checkInternalBasisFunctionsProduct(indexer, p);
        }
        else
        {
            checkInternalBasisFunctionsTrunk(indexer, p);
        }
    }

    void checkInternalBasisFunctionsProduct(const BasisFunctionIndexer& indexer, uint32_t p)
    {
        for (int elementIdx = 0; elementIdx < 4; elementIdx++)
        {
            if (elements.at(elementIdx).size() == 4)
            {
                checkInternalBasisFunctionsProductQuad(elementIdx, indexer, p);
            }
            else
            {
                checkInternalBasisFunctionsProductTri(elementIdx, indexer, p);
            }
        }
    }
    
    void checkInternalBasisFunctionsProductQuad(uint32_t elementIdx, const BasisFunctionIndexer& indexer, uint32_t p)
    {
        for (int k = 2; k <= p; k++)
        {
            for (int l = 2; l <= p; l++)
            {
                uint32_t shapeFunctionIdx = 4 + 4*(p-1) + (k-2)*(p-1) + l-2;
                BasisFunctionDescriptor desc = indexer.getBasisFunctionDescriptor(elementIdx, shapeFunctionIdx);
                uint32_t basisFunctionIdx = indexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
                EXPECT_EQ(desc, BasisFunctionDescriptor(InternalBasisFunctionDescriptor(elementIdx, k, l)));
                EXPECT_EQ(basisFunctionIdx, indexer.getBasisFunctionIndex(desc));
                EXPECT_EQ(desc, indexer.getBasisFunctionDescriptor(basisFunctionIdx));
                basisFunctionIndices.insert(basisFunctionIdx);
            }
        }
    }

    void checkInternalBasisFunctionsProductTri(uint32_t elementIdx, const BasisFunctionIndexer& indexer, int p)
    {
        for (int k = 0; k <= p-2; k++)
        {
            for (int l = 0; l <= p-2; l++)
            {
                uint32_t shapeFunctionIdx = 3 + 3*(p-1) + k*(p-1) + l;
                BasisFunctionDescriptor desc = indexer.getBasisFunctionDescriptor(elementIdx, shapeFunctionIdx);
                uint32_t basisFunctionIdx = indexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
                EXPECT_EQ(desc, BasisFunctionDescriptor(InternalBasisFunctionDescriptor(elementIdx, k, l)));
                EXPECT_EQ(basisFunctionIdx, indexer.getBasisFunctionIndex(desc));
                EXPECT_EQ(desc, indexer.getBasisFunctionDescriptor(basisFunctionIdx));
                basisFunctionIndices.insert(basisFunctionIdx);
            }
        }
    }

    void checkInternalBasisFunctionsTrunk(const BasisFunctionIndexer& indexer, int p)
    {
        for (int elementIdx = 0; elementIdx < 4; elementIdx++)
        {
            if (elements.at(elementIdx).size() == 4)
            {
                checkInternalBasisFunctionsTrunkQuad(elementIdx, indexer, p);
            }
            else
            {
                checkInternalBasisFunctionsTrunkTri(elementIdx, indexer, p);
            }
        }
    }

    void checkInternalBasisFunctionsTrunkQuad(uint32_t elementIdx, const BasisFunctionIndexer& indexer, int p)
    {
        uint32_t internalShapeFunctionIdxBase = 0;
        for (int k = 2; k <= p-2; k++)
        {
            for (int l = 2; l <= p-k; l++)
            {
                uint32_t shapeFunctionIdx = 4 + 4*(p-1) + internalShapeFunctionIdxBase + l-2;
                BasisFunctionDescriptor desc = indexer.getBasisFunctionDescriptor(elementIdx, shapeFunctionIdx);
                uint32_t basisFunctionIdx = indexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
                EXPECT_EQ(desc, BasisFunctionDescriptor(InternalBasisFunctionDescriptor(elementIdx, k, l)));
                EXPECT_EQ(basisFunctionIdx, indexer.getBasisFunctionIndex(desc));
                EXPECT_EQ(desc, indexer.getBasisFunctionDescriptor(basisFunctionIdx));
                basisFunctionIndices.insert(basisFunctionIdx);
            }
            internalShapeFunctionIdxBase += p-k-1;
        }
    }

    void checkInternalBasisFunctionsTrunkTri(uint32_t elementIdx, const BasisFunctionIndexer& indexer, int p)
    {
        uint32_t internalShapeFunctionIdxBase = 0;
        for (int k = 0; k <= p-3; k++)
        {
            for (int l = 0; l <= p-k-3; l++)
            {
                uint32_t shapeFunctionIdx = 3 + 3*(p-1) + internalShapeFunctionIdxBase + l;
                BasisFunctionDescriptor desc = indexer.getBasisFunctionDescriptor(elementIdx, shapeFunctionIdx);
                uint32_t basisFunctionIdx = indexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
                EXPECT_EQ(desc, BasisFunctionDescriptor(InternalBasisFunctionDescriptor(elementIdx, k, l)));
                EXPECT_EQ(basisFunctionIdx, indexer.getBasisFunctionIndex(desc));
                EXPECT_EQ(desc, indexer.getBasisFunctionDescriptor(basisFunctionIdx));
                basisFunctionIndices.insert(basisFunctionIdx);
            }
            internalShapeFunctionIdxBase += p-k-2;
        }
    }

    void checkContinuityOverSides(uint32_t p)
    {
        BasisFunctionIndexer indexer(mesh, p, PolynomialSpaceType_Product);
        checkContinuityOverSide2(indexer, p);
        checkContinuityOverSide3(indexer, p);
        checkContinuityOverSide5(indexer, p);
        checkContinuityOverSide8(indexer, p);
    }

    void checkContinuityOverSide2(const BasisFunctionIndexer& indexer, uint32_t p)
    {
        const uint32_t elementIdx1 = 0;
        const uint32_t elementIdx2 = 1;
        const Polynomial1D sideParameterizationX("t");
        const Polynomial1D sideParameterizationY(0);
        const AffineMap G1 = mesh.getElement(elementIdx1).getReferenceElementMap().inverse();
        const AffineMap G2 = mesh.getElement(elementIdx2).getReferenceElementMap().inverse();
        
        /* Nodal basis functions */
        EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 2), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(indexer.getShapeFunction(elementIdx2, 0), G2), sideParameterizationX, sideParameterizationY));
        EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 3), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(indexer.getShapeFunction(elementIdx2, 2), G2), sideParameterizationX, sideParameterizationY));

        /* Side basis functions */
        for (int i = 0; i < p-1; i++)
        {
            EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 4 + 2*(p-1) + i), G1), sideParameterizationX, sideParameterizationY),
                      compose(compose(indexer.getShapeFunction(elementIdx2, 3 + 2*(p-1) + i), G2), sideParameterizationX, sideParameterizationY));
        }
    }
    
    void checkContinuityOverSide3(const BasisFunctionIndexer& indexer, uint32_t p)
    {
        const uint32_t elementIdx1 = 0;
        const uint32_t elementIdx2 = 3;
        const Polynomial1D sideParameterizationX(0);
        const Polynomial1D sideParameterizationY("t-1");
        const AffineMap G1 = mesh.getElement(elementIdx1).getReferenceElementMap().inverse();
        const AffineMap G2 = mesh.getElement(elementIdx2).getReferenceElementMap().inverse();
        
        /* Nodal basis functions */
        EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 0), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(indexer.getShapeFunction(elementIdx2, 2), G2), sideParameterizationX, sideParameterizationY));
        EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 3), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(indexer.getShapeFunction(elementIdx2, 0), G2), sideParameterizationX, sideParameterizationY));

        /* Side basis functions */
        for (int i = 0; i < p-1; i++)
        {
            EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 4 + 3*(p-1) + i), G1), sideParameterizationX, sideParameterizationY),
                      compose(compose(indexer.getShapeFunction(elementIdx2, 3 + 2*(p-1) + i), G2), sideParameterizationX, sideParameterizationY));
        }
    }

    void checkContinuityOverSide5(const BasisFunctionIndexer& indexer, uint32_t p)
    {
        const uint32_t elementIdx1 = 1;
        const uint32_t elementIdx2 = 2;
        const Polynomial1D sideParameterizationX(0);
        const Polynomial1D sideParameterizationY("t");
        const AffineMap G1 = mesh.getElement(elementIdx1).getReferenceElementMap().inverse();
        const AffineMap G2 = mesh.getElement(elementIdx2).getReferenceElementMap().inverse();
        
        /* Nodal basis functions */
        EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 1), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(indexer.getShapeFunction(elementIdx2, 0), G2), sideParameterizationX, sideParameterizationY));
        EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 2), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(indexer.getShapeFunction(elementIdx2, 3), G2), sideParameterizationX, sideParameterizationY));

        /* Side basis functions */
        for (int i = 0; i < p-1; i++)
        {
            EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 3 + 1*(p-1) + i), G1), sideParameterizationX, sideParameterizationY),
                      compose(compose(indexer.getShapeFunction(elementIdx2, 4 + 3*(p-1) + i), G2), sideParameterizationX, sideParameterizationY));
        }
    }

    void checkContinuityOverSide8(const BasisFunctionIndexer& indexer, uint32_t p)
    {
        const uint32_t elementIdx1 = 2;
        const uint32_t elementIdx2 = 3;
        const Polynomial1D sideParameterizationX("-t");
        const Polynomial1D sideParameterizationY(0);
        const AffineMap G1 = mesh.getElement(elementIdx1).getReferenceElementMap().inverse();
        const AffineMap G2 = mesh.getElement(elementIdx2).getReferenceElementMap().inverse();
        
        /* Nodal basis functions */
        EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 2), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(indexer.getShapeFunction(elementIdx2, 1), G2), sideParameterizationX, sideParameterizationY));
        EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 3), G1), sideParameterizationX, sideParameterizationY),
                  compose(compose(indexer.getShapeFunction(elementIdx2, 0), G2), sideParameterizationX, sideParameterizationY));

        /* Side basis functions */
        for (int i = 0; i < p-1; i++)
        {
            EXPECT_EQ(compose(compose(indexer.getShapeFunction(elementIdx1, 4 + 2*(p-1) + i), G1), sideParameterizationX, sideParameterizationY),
                      compose(compose(indexer.getShapeFunction(elementIdx2, 3 + 0*(p-1) + i), G2), sideParameterizationX, sideParameterizationY));
        }
    }
};

TEST_F(BasisFunctionIndexerTest, ProductSpace)
{
    for (int p = 1; p <= 7; p++)
    {
        checkBasisFunctions(p, PolynomialSpaceType_Product);
    }
}

TEST_F(BasisFunctionIndexerTest, TrunkSpace)
{
    for (int p = 1; p <= 7; p++)
    {
        checkBasisFunctions(p, PolynomialSpaceType_Trunk);
    }
}

TEST_F(BasisFunctionIndexerTest, ContinuityOverSides)
{
    for (int p = 1; p <= 7; p++)
    {
        checkContinuityOverSides(p);
    }
}
} // namespace fem::ut
