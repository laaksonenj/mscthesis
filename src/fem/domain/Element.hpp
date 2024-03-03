#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "fem/domain/Node.hpp"
#include "fem/domain/Side.hpp"
#include "fem/math/AffineMap.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
enum ElementType
{
    ElementType_Triangle,
    ElementType_Parallelogram
};

class Element
{
public:
    virtual ~Element() = default;

    virtual ElementType getElementType() const = 0;
    virtual uint32_t getNumOfNodes() const = 0;
    virtual Node getNode(uint32_t nodeIdx) const = 0;
    virtual AffineMap getReferenceElementMap() const = 0;
    uint32_t getNumOfSides() const { return getNumOfNodes(); }
    Side getSide(uint32_t sideIdx) const;
    virtual std::vector<std::unique_ptr<Element>> subdivide(const Vector2mpq& x) const = 0;
};

bool areIntersecting(const Element& element1, const Element& element2);
bool isIntersectionOneNode(const Element& element1, const Element& element2);
bool isIntersectionOneSide(const Element& element1, const Element& element2);
bool isPointInsideElement(const Vector2mpq& point, const Element& element);
bool operator==(const Element& lhs, const Element& rhs);
bool operator!=(const Element& lhs, const Element& rhs);
} // namespace fem

#include "fem/domain/Parallelogram.hpp"
#include "fem/domain/Triangle.hpp"
