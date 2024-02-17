#include <gtest/gtest.h>

#include <sstream>

#include "fem/domain/Mesh.hpp"

namespace fem::ut
{
namespace
{
std::vector<Node> getFixtureNodes()
{
    return {{"-1", "1"},
            {"-1", "0"},
            {"0", "1"},
            {"0", "0"},
            {"-1", "-1"},
            {"0", "-1"},
            {"1/2", "-1"},
            {"1", "-1/2"},
            {"1", "0"},
            {"1", "1/2"},
            {"1/2", "0"},
            {"1/2", "-1/2"}};
}

std::vector<std::vector<Mesh::NodeIndex>> getFixtureElements()
{
    return {{0,1,2},
            {1,3,2},
            {1,4,5,3},
            {11,10,3},
            {3,5,11},
            {6,11,5},
            {6,7,8,11},
            {9,10,11,8}};
}
} // namespace

class MeshTestFixture : public testing::Test
{
protected:
    static inline const std::vector<Node> nodes{getFixtureNodes()};
    static inline const std::vector<std::vector<Mesh::NodeIndex>> elements{getFixtureElements()};
    static inline const Mesh mesh{nodes, elements};
};

TEST_F(MeshTestFixture, GetNumOfNodes)
{
    EXPECT_EQ(mesh.getNumOfNodes(), nodes.size());
}

TEST_F(MeshTestFixture, GetNumOfSides)
{
    EXPECT_EQ(mesh.getNumOfSides(), 19);
}

TEST_F(MeshTestFixture, GetNumOfElements)
{
    EXPECT_EQ(mesh.getNumOfElements(), elements.size());
}

TEST_F(MeshTestFixture, GetGlobalNodeIndex)
{
    for (int elementIdx = 0; elementIdx < elements.size(); elementIdx++)
    {
        const auto& globalVertexIndices = elements[elementIdx];
        for (int localVertexIdx = 0; localVertexIdx < globalVertexIndices.size(); localVertexIdx++)
        {
            EXPECT_EQ(mesh.getGlobalNodeIndex(elementIdx, localVertexIdx), globalVertexIndices[localVertexIdx]);
        }
    }
}

TEST_F(MeshTestFixture, GetGlobalSideIndex)
{
    EXPECT_EQ(mesh.getGlobalSideIndex(1, 0), 3);
    EXPECT_EQ(mesh.getGlobalSideIndex(1, 1), 4);
    EXPECT_EQ(mesh.getGlobalSideIndex(1, 2), 1);

    EXPECT_EQ(mesh.getGlobalSideIndex(4, 0), 7);
    EXPECT_EQ(mesh.getGlobalSideIndex(4, 1), 11);
    EXPECT_EQ(mesh.getGlobalSideIndex(4, 2), 10);

    EXPECT_EQ(mesh.getGlobalSideIndex(7, 0), 17);
    EXPECT_EQ(mesh.getGlobalSideIndex(7, 1), 8);
    EXPECT_EQ(mesh.getGlobalSideIndex(7, 2), 16);
    EXPECT_EQ(mesh.getGlobalSideIndex(7, 3), 18);
}

TEST_F(MeshTestFixture, GetElement)
{
    for (int elementIdx = 0; elementIdx < elements.size(); elementIdx++)
    {
        const Element& element = mesh.getElement(elementIdx);
        const auto& idx = elements[elementIdx];
        if (idx.size() == 3)
        {
            EXPECT_EQ(element, Triangle(nodes.at(idx[0]), nodes.at(idx[1]), nodes.at(idx[2])));
        }
        else if (idx.size() == 4)
        {
            EXPECT_EQ(element, Parallelogram(nodes.at(idx[0]), nodes.at(idx[1]), nodes.at(idx[2]), nodes.at(idx[3])));
        }
        else
        {
            FAIL();
        }
    }
}

TEST_F(MeshTestFixture, GetIndexOfAdjacentElement)
{
    EXPECT_FALSE(mesh.getIndexOfAdjacentElement(0, 0).has_value());
    ASSERT_TRUE(mesh.getIndexOfAdjacentElement(0, 1).has_value());
    EXPECT_EQ(mesh.getIndexOfAdjacentElement(0, 1).value(), 1);
    EXPECT_FALSE(mesh.getIndexOfAdjacentElement(0, 2).has_value());

    ASSERT_TRUE(mesh.getIndexOfAdjacentElement(4, 0).has_value());
    EXPECT_EQ(mesh.getIndexOfAdjacentElement(4, 0).value(), 2);
    ASSERT_TRUE(mesh.getIndexOfAdjacentElement(4, 1).has_value());
    EXPECT_EQ(mesh.getIndexOfAdjacentElement(4, 1).value(), 5);
    ASSERT_TRUE(mesh.getIndexOfAdjacentElement(4, 2).has_value());
    EXPECT_EQ(mesh.getIndexOfAdjacentElement(4, 2).value(), 3);

    EXPECT_FALSE(mesh.getIndexOfAdjacentElement(6, 0).has_value());
    EXPECT_FALSE(mesh.getIndexOfAdjacentElement(6, 1).has_value());
    ASSERT_TRUE(mesh.getIndexOfAdjacentElement(6, 2).has_value());
    EXPECT_EQ(mesh.getIndexOfAdjacentElement(6, 2).value(), 7);
    ASSERT_TRUE(mesh.getIndexOfAdjacentElement(6, 3).has_value());
    EXPECT_EQ(mesh.getIndexOfAdjacentElement(6, 3).value(), 5);
}

TEST_F(MeshTestFixture, GetIndexOfElementContainingPoint)
{
    EXPECT_EQ(getIndexOfElementContainingPoint(mesh, Vector2mpq(mpq_class("3/4"), 0)), 7);
    EXPECT_EQ(getIndexOfElementContainingPoint(mesh, Vector2mpq(0, mpq_class("-1/2"))), 2);
    EXPECT_EQ(getIndexOfElementContainingPoint(mesh, Vector2mpq(mpq_class("1/2"), mpq_class("-1/2"))), 3);
}

TEST_F(MeshTestFixture, GetMeshBoundary)
{
    const std::vector<std::pair<Mesh::ElementIndex, Mesh::SideIndex>> expected = {
        {0, 0}, {0, 2}, {1, 1}, {2, 0}, {2, 1}, {3, 1}, {5, 2}, {6, 0}, {6, 1}, {7, 0}, {7, 3}
    };
    EXPECT_EQ(getMeshBoundary(mesh), expected);
}

TEST(MeshTest, ParseMeshFile)
{
    std::stringstream ss;
    ss << "n " << mpq_class(-1) << " " << mpq_class(-1) << '\n';
    ss << "n " << mpq_class("1/3") << " " << mpq_class(-1) << '\n';
    ss << "n " << mpq_class("1/3") << " " << mpq_class(10) << '\n';
    ss << "n " << mpq_class(-1) << " " << mpq_class(10) << '\n';
    ss << "n " << mpq_class("-1/3") << " " << mpq_class("101/10") << '\n';
    ss << "e " << 0 << " " << 1 << " " << 2 << " " << 3 << '\n';
    ss << "e " << 2 << " " << 4 << " " << 3 << '\n';
    const auto nodesAndElements = parseMeshFile(ss);
    EXPECT_EQ(nodesAndElements.first, (std::vector<Node>{{-1, -1}, {"1/3", -1}, {"1/3", 10}, {-1, 10}, {"-1/3", "101/10"}}));
    EXPECT_EQ(nodesAndElements.second, (std::vector<std::vector<Mesh::NodeIndex>>{{0,1,2,3}, {2,4,3}}));
}

TEST(MeshTest, NonConformingMeshAborts)
{
    std::vector<Node> nodes;
    std::vector<std::vector<Mesh::NodeIndex>> elements;

    /* Overlapping elements */
    nodes = getFixtureNodes();
    elements = getFixtureElements();
    int s = nodes.size();
    nodes.push_back({"-3/4", "3/4"});
    elements.push_back({2,0,s});
    EXPECT_DEATH(Mesh(nodes, elements), ".*");

    /* An edge is only partially shared */
    nodes = getFixtureNodes();
    elements = getFixtureElements();
    s = nodes.size();
    nodes.push_back({"-1/2", "0"});
    nodes.push_back({"-1/2", "-1"});
    elements.erase(elements.begin() + 2);
    elements.push_back({s+1,s,1,4});
    elements.push_back({s,s+1,5,3});
    EXPECT_DEATH(Mesh(nodes, elements), ".*");
}
} // namespace fem::ut
