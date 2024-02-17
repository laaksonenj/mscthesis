#pragma once

#include <array>
#include <cassert>

#include "fem/domain/Element.hpp"

namespace fem
{
class Triangle : public Element
{
public:
    explicit Triangle(const Node& n1, const Node& n2, const Node& n3);

    ElementType getElementType() const override { return ElementType_Triangle; };
    uint32_t getNumOfNodes() const override { return 3; };
    Node getNode(uint32_t nodeIdx) const override { assert(nodeIdx < 3); return m_nodes[nodeIdx]; }
    AffineMap getReferenceElementMap() const override;

private:
    std::array<Node, 3> m_nodes;
};
} // namespace fem
