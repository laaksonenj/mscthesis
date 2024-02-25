#pragma once

#include <cassert>
#include <cstdint>
#include <istream>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "fem/domain/Element.hpp"
#include "fem/domain/Node.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
class Mesh
{
public:
    using NodeIndex = int;
    using SideIndex = int;
    using ElementIndex = int;

public:
    explicit Mesh(const std::vector<Node>& nodes, const std::vector<std::vector<NodeIndex>>& elements);

    uint32_t getNumOfNodes() const { return m_numOfNodes; }
    uint32_t getNumOfSides() const { return m_numOfSides; }
    uint32_t getNumOfElements() const { return m_elements.size(); }

    NodeIndex getGlobalNodeIndex(ElementIndex elementIdx, NodeIndex localNodeIdx) const;
    SideIndex getGlobalSideIndex(ElementIndex elementIdx, SideIndex localSideIdx) const;

    const Element& getElement(ElementIndex elementIdx) const { assert(elementIdx < getNumOfElements()); return *(m_elements[elementIdx]); }
    std::optional<ElementIndex> getIndexOfAdjacentElement(ElementIndex elementIdx, SideIndex localSideIdx) const;

    bool containsTriangle() const { return m_containsTriangle; }
    bool containsQuadrilateral() const { return m_containsQuadrilateral; }

private:
    void createElements(const std::vector<Node>& nodes, const std::vector<std::vector<NodeIndex>>& elements);
    void assignGlobalSideIndices();
    void updateAdjacentElements();
    void validateElements();

private:
    std::vector<std::unique_ptr<Element>> m_elements;
    std::vector<std::vector<NodeIndex>> m_globalNodeIndices;
    std::vector<std::vector<SideIndex>> m_globalSideIndices;
    std::vector<std::vector<ElementIndex>> m_adjacentElements;
    uint32_t m_numOfNodes;
    uint32_t m_numOfSides;
    bool m_containsTriangle;
    bool m_containsQuadrilateral;
};

Mesh::ElementIndex getIndexOfElementContainingPoint(const Mesh& mesh, const Vector2mpq& point);
const Element& getElementContainingPoint(const Mesh& mesh, const Vector2mpq& point);
std::vector<std::pair<Mesh::ElementIndex, Mesh::SideIndex>> getMeshBoundary(const Mesh& mesh);
Mesh createMeshFromFile(const std::string& filename);
Mesh createMeshFromFile(std::istream& input);
std::pair<std::vector<Node>, std::vector<std::vector<Mesh::NodeIndex>>> parseMeshFile(std::istream& input);
mpq_class calculateMeshArea(const Mesh& mesh);
} // namespace fem
