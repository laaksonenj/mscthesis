#include "fem/domain/Mesh.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <utility>

namespace fem
{
Mesh::Mesh(const std::vector<Node>& nodes, const std::vector<std::vector<NodeIndex>>& elements)
    : m_globalNodeIndices(elements)
    , m_numOfNodes(nodes.size())
    , m_containsTriangle(false)
    , m_containsQuadrilateral(false)
{
    createElements(nodes, elements);
    assignGlobalSideIndices();
    updateAdjacentElements();
    validateElements();
}

Mesh::NodeIndex Mesh::getGlobalNodeIndex(ElementIndex elementIdx, NodeIndex localNodeIdx) const
{
    assert(elementIdx < m_globalNodeIndices.size());
    const std::vector<NodeIndex>& indices = m_globalNodeIndices[elementIdx];
    assert(localNodeIdx < indices.size());
    return indices[localNodeIdx];
}

Mesh::SideIndex Mesh::getGlobalSideIndex(ElementIndex elementIdx, SideIndex localSideIdx) const
{
    assert(elementIdx < m_globalSideIndices.size());
    const std::vector<SideIndex>& indices = m_globalSideIndices[elementIdx];
    assert(localSideIdx < indices.size());
    return indices[localSideIdx];
}

std::optional<Mesh::ElementIndex> Mesh::getIndexOfAdjacentElement(ElementIndex elementIdx, SideIndex localSideIdx) const
{
    assert(elementIdx < m_adjacentElements.size());
    const std::vector<ElementIndex>& indices = m_adjacentElements[elementIdx];
    assert(localSideIdx < indices.size());
    if (indices[localSideIdx] != -1)
    {
        return indices[localSideIdx];
    }
    else
    {
        return {};
    }
}

void Mesh::createElements(const std::vector<Node>& nodes, const std::vector<std::vector<NodeIndex>>& elements)
{
    for (const auto& idxs : elements)
    {
        if (idxs.size() == 3)
        {
            m_elements.push_back(std::make_unique<Triangle>(nodes.at(idxs[0]), nodes.at(idxs[1]), nodes.at(idxs[2])));
            m_containsTriangle = true;
        }
        else if (idxs.size() == 4)
        {
            m_elements.push_back(std::make_unique<Parallelogram>(nodes.at(idxs[0]), nodes.at(idxs[1]), nodes.at(idxs[2]), nodes.at(idxs[3])));
            m_containsQuadrilateral = true;
        }
        else
        {
            assert(false && "Invalid element");
        }
    }
}

void Mesh::assignGlobalSideIndices()
{
    std::map<std::pair<NodeIndex, NodeIndex>, SideIndex> sides;
    SideIndex sideIdx = 0;
    for (const auto& nodeIdxs : m_globalNodeIndices)
    {
        std::vector<SideIndex> sideIdxs;
        for (int i = 0; i < nodeIdxs.size(); i++)
        {
            const NodeIndex a = nodeIdxs[i];
            const NodeIndex b = nodeIdxs[(i + 1) % nodeIdxs.size()];
            const auto key = std::make_pair(std::min(a, b), std::max(a, b));
            if (!sides.contains(key))
            {
                sides.emplace(key, sideIdx);
                sideIdx++;
            }
            sideIdxs.push_back(sides.at(key));
        }
        m_globalSideIndices.push_back(sideIdxs);
    }
    m_numOfSides = sideIdx;
}

void Mesh::updateAdjacentElements()
{
    std::map<SideIndex, ElementIndex> elementAdjacentToSideMap;
    for (int elementIdx = 0; elementIdx < m_globalSideIndices.size(); elementIdx++)
    {
        std::vector<ElementIndex> elementIdxs;
        for (const auto sideIdx : m_globalSideIndices[elementIdx])
        {
            if (elementAdjacentToSideMap.contains(sideIdx))
            {
                const ElementIndex adjacentElementIdx = elementAdjacentToSideMap.at(sideIdx);
                assert(adjacentElementIdx < m_adjacentElements.size());
                assert(adjacentElementIdx < m_globalSideIndices.size());
                auto it = std::find(m_globalSideIndices[adjacentElementIdx].begin(), m_globalSideIndices[adjacentElementIdx].end(), sideIdx);
                assert(m_adjacentElements[adjacentElementIdx][std::distance(m_globalSideIndices[adjacentElementIdx].begin(), it)] == -1);
                m_adjacentElements[adjacentElementIdx][std::distance(m_globalSideIndices[adjacentElementIdx].begin(), it)] = elementIdx;
                elementIdxs.push_back(adjacentElementIdx);
            }
            else
            {
                elementAdjacentToSideMap.emplace(sideIdx, elementIdx);
                elementIdxs.push_back(-1);
            }
        }
        m_adjacentElements.push_back(elementIdxs);
    }
}

void Mesh::validateElements()
{
    for (int i = 0; i < getNumOfElements(); i++)
    {
        const Element& element1 = getElement(i);
        for (int j = i + 1; j < getNumOfElements(); j++)
        {
            const Element& element2 = getElement(j);
            assert(!areIntersecting(element1, element2) || isIntersectionOneNode(element1, element2) || isIntersectionOneSide(element1, element2));
        }
    }
}

Mesh::ElementIndex getIndexOfElementContainingPoint(const Mesh& mesh, const Vector2mpq& point)
{
    for (int elementIdx = 0; elementIdx < mesh.getNumOfElements(); elementIdx++)
    {
        if (isPointInsideElement(point, mesh.getElement(elementIdx)))
        {
            return elementIdx;
        }
    }
    assert(false && "Point is outside the mesh");
    return -1;
}

const Element& getElementContainingPoint(const Mesh& mesh, const Vector2mpq& point)
{
    return mesh.getElement(getIndexOfElementContainingPoint(mesh, point));
}

std::vector<std::pair<Mesh::ElementIndex, Mesh::SideIndex>> getMeshBoundary(const Mesh& mesh)
{
    std::vector<std::pair<Mesh::ElementIndex, Mesh::SideIndex>> res;
    for (int elementIdx = 0; elementIdx < mesh.getNumOfElements(); elementIdx++)
    {
        const auto& element = mesh.getElement(elementIdx);
        for (int sideIdx = 0; sideIdx < element.getNumOfSides(); sideIdx++)
        {
            if (!mesh.getIndexOfAdjacentElement(elementIdx, sideIdx).has_value())
            {
                res.push_back(std::make_pair(elementIdx, sideIdx));
            }
        }
    }
    return res;
}

Mesh createMeshFromFile(const std::string& filename)
{
    std::ifstream ifs(filename);
    if (!ifs.is_open())
    {
        std::cout << "Could not open the file " << filename << std::endl;
        assert(false);
    }
    return createMeshFromFile(ifs);
}

Mesh createMeshFromFile(std::istream& input)
{
    const auto nodesAndElements = parseMeshFile(input);
    return Mesh(nodesAndElements.first, nodesAndElements.second);
}

std::pair<std::vector<Node>, std::vector<std::vector<Mesh::NodeIndex>>> parseMeshFile(std::istream& input)
{
    std::vector<Node> nodes;
    std::vector<std::vector<Mesh::NodeIndex>> elements;
    for (std::string line; std::getline(input, line);)
    {
        std::stringstream ss(line);
        std::string token;
        std::getline(ss, token, ' ');
        if (token == "n")
        {
            std::string x, y;
            std::getline(ss, x, ' ');
            std::getline(ss, y, ' ');
            if (x.empty() || y.empty())
            {
                std::cout << "Invalid mesh file line: " << line << std::endl;
                assert(false);
            }
            nodes.push_back(Node(mpq_class(x), mpq_class(y)));
        }
        else if (token == "e")
        {
            std::vector<int> e;
            int i;
            while (ss >> i)
            {
                e.push_back(i);
            }
            elements.push_back(e);
        }
        else
        {
            std::cout << "Invalid mesh file line: " << line << std::endl;
            assert(false);
        }
    }
    return std::make_pair(std::move(nodes), std::move(elements));
}

mpq_class calculateMeshArea(const Mesh& mesh)
{
    mpq_class res = 0;
    for (int elementIdx = 0; elementIdx < mesh.getNumOfElements(); elementIdx++)
    {
        const Element& element = mesh.getElement(elementIdx);
        const mpq_class detA = element.getReferenceElementMap().A.determinant();
        if (element.getElementType() == ElementType_Parallelogram)
        {
            res += 4 * detA;
        }
        else
        {
            res += detA / 2;
        }
    }
    return res;
}
} // namespace fem
