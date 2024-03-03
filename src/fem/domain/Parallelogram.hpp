#pragma once

#include <array>
#include <cassert>

#include "fem/domain/Element.hpp"

namespace fem
{
class Parallelogram : public Element
{
public:
    explicit Parallelogram(const Node& n1, const Node& n2, const Node& n3, const Node& n4);

    ElementType getElementType() const override { return ElementType_Parallelogram; }
    uint32_t getNumOfNodes() const override { return 4; }
    Node getNode(uint32_t nodeIdx) const override { assert(nodeIdx < 4); return m_nodes[nodeIdx]; }
    AffineMap getReferenceElementMap() const override;
    std::vector<std::unique_ptr<Element>> subdivide(const Vector2mpq& x) const override;

private:
    std::array<Node, 4> m_nodes;
};
} // namespace fem
