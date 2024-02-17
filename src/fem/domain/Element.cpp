#include "fem/domain/Element.hpp"

#include <algorithm>
#include <cassert>
#include <vector>

namespace fem
{
Side Element::getSide(uint32_t sideIdx) const
{
    assert(sideIdx < getNumOfSides());
    const Node a = getNode(sideIdx);
    const Node b = getNode((sideIdx + 1) % getNumOfNodes());
    return Side(a, b);
}

bool operator==(const Element& lhs, const Element& rhs)
{
    if (lhs.getNumOfNodes() != rhs.getNumOfNodes())
    {
        return false;
    }
    for (int i = 0; i < lhs.getNumOfNodes(); i++)
    {
        if (lhs.getNode(i) != rhs.getNode(i))
        {
            return false;
        }
    }
    return true;
}

bool operator!=(const Element& lhs, const Element& rhs)
{
    return !(lhs == rhs);
}

bool isPointInsideElement(const Vector2mpq& point, const Element& element)
{
    const AffineMap invMap = element.getReferenceElementMap().inverse();
    const Vector2mpq vertexInLocalCoords = invMap(point);
    const mpq_class& x = vertexInLocalCoords(0);
    const mpq_class& y = vertexInLocalCoords(1);
    if (element.getElementType() == ElementType_Triangle)
    {
        return x >= 0 && y >= 0 && x + y <= 1;
    }
    else if (element.getElementType() == ElementType_Parallelogram)
    {
        return x >= -1 && x <= 1 && y >= -1 && y <= 1;
    }
    else
    {
        assert(false);
        return false;
    }
}

bool areIntersecting(const Element& element1, const Element& element2)
{
    for (int i = 0; i < element1.getNumOfNodes(); i++)
    {
        if (isPointInsideElement(element1.getNode(i), element2))
        {
            return true;
        }
    }
    for (int i = 0; i < element2.getNumOfNodes(); i++)
    {
        if (isPointInsideElement(element2.getNode(i), element1))
        {
            return true;
        }
    }
    for (int i = 0; i < element1.getNumOfSides(); i++)
    {
        for (int j = 0; j < element2.getNumOfSides(); j++)
        {
            if (areIntersecting(element1.getSide(i), element2.getSide(j)))
            {
                return true;
            }
        }
    }
    return false;
}

bool isIntersectionOneNode(const Element& element1, const Element& element2)
{
    bool hasOneCommonNode = false;
    int nodeIndex1, nodeIndex2;
    for (int i = 0; i < element1.getNumOfNodes(); i++)
    {
        const Node n1 = element1.getNode(i);
        for (int j = 0; j < element2.getNumOfNodes(); j++)
        {
            const Node n2 = element2.getNode(j);
            if (n1 == n2)
            {
                if (hasOneCommonNode)
                {
                    return false;
                }
                hasOneCommonNode = true;
                nodeIndex1 = i;
                nodeIndex2 = j;
            }
        }
    }
    if (!hasOneCommonNode)
    {
        return false;
    }
    const int i1 = (nodeIndex1 + element1.getNumOfNodes()- 1) % element1.getNumOfNodes();
    const int j1 = (nodeIndex1 + 1) % element1.getNumOfNodes();
    const int i2 = (nodeIndex2 + element2.getNumOfNodes() - 1) % element2.getNumOfNodes();
    const int j2 = (nodeIndex2 + 1) % element2.getNumOfNodes();
    if (isPointInsideElement(element1.getNode(i1), element2) ||
        isPointInsideElement(element1.getNode(j1), element2) ||
        isPointInsideElement(element2.getNode(i2), element1) ||
        isPointInsideElement(element2.getNode(j2), element1))
    {
        return false;
    }
    const Side side1 = Side(element1.getNode(i1), element1.getNode(j1));
    const Side side2 = Side(element2.getNode(i2), element2.getNode(j2));
    if (areIntersecting(side1, element2.getSide(i2)) ||
        areIntersecting(side1, element2.getSide(nodeIndex2)) ||
        areIntersecting(side2, element1.getSide(i1)) ||
        areIntersecting(side2, element1.getSide(nodeIndex1)))
    {
        return false;
    }
    return true;
}

bool isIntersectionOneSide(const Element& element1, const Element& element2)
{
    for (int i = 0; i < element1.getNumOfSides(); i++)
    {
        const Side side1 = element1.getSide(i);
        for (int j = 0; j < element2.getNumOfSides(); j++)
        {
            const Side side2 = element2.getSide(j);
            if (side1 == side2)
            {
                if (element1.getNode(i) == element2.getNode(j))
                {
                    return false;
                }
                else
                {
                    return true;
                }
            }
        }
    }
    return false;
}
} // namespace fem
