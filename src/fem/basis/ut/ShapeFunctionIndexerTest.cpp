#include <gtest/gtest.h>

#include <unordered_set>

#include "fem/basis/ShapeFunctionIndexer.hpp"

namespace fem::ut
{
TEST(ShapeFunctionIndexerTest, TriangleTrunkSpace)
{
    const ElementType et = ElementType_Triangle;
    for (int p = 1; p <= 7; p++)
    {
        ShapeFunctionIndexer indexer(p, PolynomialSpaceType_Trunk);
        EXPECT_EQ(indexer.getNumOfNodalShapeFunctions(et), 3);
        EXPECT_EQ(indexer.getNumOfSideShapeFunctions(et), 3*(p-1));
        EXPECT_EQ(indexer.getNumOfInternalShapeFunctions(et), (p-1)*(p-2)/2);
        EXPECT_EQ(indexer.getNumOfShapeFunctions(et), 3 + 3*(p-1) + (p-1)*(p-2)/2);
        
        std::unordered_set<ShapeFunctionDescriptor> descs;
        for (int i = 0; i < indexer.getNumOfShapeFunctions(et); i++)
        {
            descs.insert(indexer.getShapeFunctionDescriptor(et, i));
        }
        EXPECT_EQ(descs.size(), indexer.getNumOfShapeFunctions(et));
        for (int nodeIdx = 0; nodeIdx < 3; nodeIdx++)
        {
            EXPECT_TRUE(descs.contains(NodalShapeFunctionDescriptor(nodeIdx)));
        }
        for (int sideIdx = 0; sideIdx < 3; sideIdx++)
        {
            for (int k = 2; k <= p; k++)
            {
                EXPECT_TRUE(descs.contains(SideShapeFunctionDescriptor(sideIdx, k)));
            }
        }
        uint32_t internalShapeFunctionIndex = 0;
        for (int k = 0; k <= p-3; k++)
        {
            for (int l = 0; l <= p-3-k; l++)
            {
                const InternalShapeFunctionDescriptor desc(k, l);
                EXPECT_TRUE(descs.contains(desc));
                EXPECT_EQ(indexer.getInternalShapeFunctionIndex(et, desc), internalShapeFunctionIndex);
                EXPECT_EQ(indexer.getInternalShapeFunctionDescriptor(et, internalShapeFunctionIndex), desc);
                internalShapeFunctionIndex++;
            }
        }
    }
}

TEST(ShapeFunctionIndexerTest, TriangleProductSpace)
{
    const ElementType et = ElementType_Triangle;
    for (int p = 1; p <= 7; p++)
    {
        ShapeFunctionIndexer indexer(p, PolynomialSpaceType_Product);
        EXPECT_EQ(indexer.getNumOfNodalShapeFunctions(et), 3);
        EXPECT_EQ(indexer.getNumOfSideShapeFunctions(et), 3*(p-1));
        EXPECT_EQ(indexer.getNumOfInternalShapeFunctions(et), (p-1)*(p-1));
        EXPECT_EQ(indexer.getNumOfShapeFunctions(et), 3 + 3*(p-1) + (p-1)*(p-1));
        
        std::unordered_set<ShapeFunctionDescriptor> descs;
        for (int i = 0; i < indexer.getNumOfShapeFunctions(et); i++)
        {
            descs.insert(indexer.getShapeFunctionDescriptor(et, i));
        }
        EXPECT_EQ(descs.size(), indexer.getNumOfShapeFunctions(et));
        for (int nodeIdx = 0; nodeIdx < 3; nodeIdx++)
        {
            EXPECT_TRUE(descs.contains(NodalShapeFunctionDescriptor(nodeIdx)));
        }
        for (int sideIdx = 0; sideIdx < 3; sideIdx++)
        {
            for (int k = 2; k <= p; k++)
            {
                EXPECT_TRUE(descs.contains(SideShapeFunctionDescriptor(sideIdx, k)));
            }
        }
        uint32_t internalShapeFunctionIndex = 0;
        for (int k = 0; k <= p-2; k++)
        {
            for (int l = 0; l <= p-2; l++)
            {
                const InternalShapeFunctionDescriptor desc(k, l);
                EXPECT_TRUE(descs.contains(desc));
                EXPECT_EQ(indexer.getInternalShapeFunctionIndex(et, desc), internalShapeFunctionIndex);
                EXPECT_EQ(indexer.getInternalShapeFunctionDescriptor(et, internalShapeFunctionIndex), desc);
                internalShapeFunctionIndex++;
            }
        }
    }
}

TEST(ShapeFunctionIndexerTest, QuadrilateralTrunkSpace)
{
    const ElementType et = ElementType_Parallelogram;
    for (int p = 1; p <= 7; p++)
    {
        ShapeFunctionIndexer indexer(p, PolynomialSpaceType_Trunk);
        EXPECT_EQ(indexer.getNumOfNodalShapeFunctions(et), 4);
        EXPECT_EQ(indexer.getNumOfSideShapeFunctions(et), 4*(p-1));
        const uint32_t numOfInternalShapeFunctions = p >= 4 ? (p-2)*(p-3)/2 : 0;
        EXPECT_EQ(indexer.getNumOfInternalShapeFunctions(et), numOfInternalShapeFunctions);
        EXPECT_EQ(indexer.getNumOfShapeFunctions(et), 4 + 4*(p-1) + numOfInternalShapeFunctions);
        
        std::unordered_set<ShapeFunctionDescriptor> descs;
        for (int i = 0; i < indexer.getNumOfShapeFunctions(et); i++)
        {
            descs.insert(indexer.getShapeFunctionDescriptor(et, i));
        }
        EXPECT_EQ(descs.size(), indexer.getNumOfShapeFunctions(et));
        for (int nodeIdx = 0; nodeIdx < 4; nodeIdx++)
        {
            EXPECT_TRUE(descs.contains(NodalShapeFunctionDescriptor(nodeIdx)));
        }
        for (int sideIdx = 0; sideIdx < 4; sideIdx++)
        {
            for (int k = 2; k <= p; k++)
            {
                EXPECT_TRUE(descs.contains(SideShapeFunctionDescriptor(sideIdx, k)));
            }
        }
        uint32_t internalShapeFunctionIndex = 0;
        for (int k = 2; k <= p-2; k++)
        {
            for (int l = 2; l <= p-k; l++)
            {
                const InternalShapeFunctionDescriptor desc(k, l);
                EXPECT_TRUE(descs.contains(desc));
                EXPECT_EQ(indexer.getInternalShapeFunctionIndex(et, desc), internalShapeFunctionIndex);
                EXPECT_EQ(indexer.getInternalShapeFunctionDescriptor(et, internalShapeFunctionIndex), desc);
                internalShapeFunctionIndex++;
            }
        }
    }
}

TEST(ShapeFunctionIndexerTest, QuadrilateralProductSpace)
{
    const ElementType et = ElementType_Parallelogram;
    for (int p = 1; p <= 7; p++)
    {
        ShapeFunctionIndexer indexer(p, PolynomialSpaceType_Product);
        EXPECT_EQ(indexer.getNumOfNodalShapeFunctions(et), 4);
        EXPECT_EQ(indexer.getNumOfSideShapeFunctions(et), 4*(p-1));
        EXPECT_EQ(indexer.getNumOfInternalShapeFunctions(et), (p-1)*(p-1));
        EXPECT_EQ(indexer.getNumOfShapeFunctions(et), 4 + 4*(p-1) + (p-1)*(p-1));
        
        std::unordered_set<ShapeFunctionDescriptor> descs;
        for (int i = 0; i < indexer.getNumOfShapeFunctions(et); i++)
        {
            descs.insert(indexer.getShapeFunctionDescriptor(et, i));
        }
        EXPECT_EQ(descs.size(), indexer.getNumOfShapeFunctions(et));
        for (int nodeIdx = 0; nodeIdx < 4; nodeIdx++)
        {
            EXPECT_TRUE(descs.contains(NodalShapeFunctionDescriptor(nodeIdx)));
        }
        for (int sideIdx = 0; sideIdx < 4; sideIdx++)
        {
            for (int k = 2; k <= p; k++)
            {
                EXPECT_TRUE(descs.contains(SideShapeFunctionDescriptor(sideIdx, k)));
            }
        }
        uint32_t internalShapeFunctionIndex = 0;
        for (int k = 2; k <= p; k++)
        {
            for (int l = 2; l <= p; l++)
            {
                const InternalShapeFunctionDescriptor desc(k, l);
                EXPECT_TRUE(descs.contains(desc));
                EXPECT_EQ(indexer.getInternalShapeFunctionIndex(et, desc), internalShapeFunctionIndex);
                EXPECT_EQ(indexer.getInternalShapeFunctionDescriptor(et, internalShapeFunctionIndex), desc);
                internalShapeFunctionIndex++;
            }
        }
    }
}
} // namespace fem::ut
