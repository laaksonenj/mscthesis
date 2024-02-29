#include "fem/basis/ShapeFunctionIndexer.hpp"

#include <cassert>

namespace fem
{
namespace
{
uint32_t calculateNumOfInternalShapeFunctionsQuadrilateral(uint32_t p, PolynomialSpaceType polynomialSpaceType)
{
    if (polynomialSpaceType == PolynomialSpaceType_Product)
    {
        return (p-1)*(p-1);
    }
    else
    {
        if (p >= 4)
        {
            return (p-2)*(p-3)/2;
        }
        else
        {
            return 0;
        }
    }
}

uint32_t calculateNumOfInternalShapeFunctionsTriangle(uint32_t p, PolynomialSpaceType polynomialSpaceType)
{
    if (polynomialSpaceType == PolynomialSpaceType_Product)
    {
        return (p-1)*(p-1);
    }
    else
    {
        return (p-1)*(p-2)/2;
    }
}
} // namespace

ShapeFunctionIndexer::ShapeFunctionIndexer(uint32_t p, PolynomialSpaceType polynomialSpaceType)
    : m_p(p)
    , m_polynomialSpaceType(polynomialSpaceType)
    , m_numOfNodalShapeFunctionsQuadrilateral(4)
    , m_numOfSideShapeFunctionsQuadrilateral(4*(p-1))
    , m_numOfInternalShapeFunctionsQuadrilateral(calculateNumOfInternalShapeFunctionsQuadrilateral(p, polynomialSpaceType))
    , m_numOfShapeFunctionsQuadrilateral(m_numOfNodalShapeFunctionsQuadrilateral + m_numOfSideShapeFunctionsQuadrilateral + m_numOfInternalShapeFunctionsQuadrilateral)
    , m_numOfNodalShapeFunctionsTriangle(3)
    , m_numOfSideShapeFunctionsTriangle(3*(p-1))
    , m_numOfInternalShapeFunctionsTriangle(calculateNumOfInternalShapeFunctionsTriangle(p, polynomialSpaceType))
    , m_numOfShapeFunctionsTriangle(m_numOfNodalShapeFunctionsTriangle + m_numOfSideShapeFunctionsTriangle + m_numOfInternalShapeFunctionsTriangle)
{
    assert(p > 0);
}

uint32_t ShapeFunctionIndexer::getNumOfNodalShapeFunctions(ElementType elementType) const
{
    if (elementType == ElementType_Parallelogram)
    {
        return m_numOfNodalShapeFunctionsQuadrilateral;
    }
    else
    {
        return m_numOfNodalShapeFunctionsTriangle;
    }
}

uint32_t ShapeFunctionIndexer::getNumOfSideShapeFunctions(ElementType elementType) const
{
    if (elementType == ElementType_Parallelogram)
    {
        return m_numOfSideShapeFunctionsQuadrilateral;
    }
    else
    {
        return m_numOfSideShapeFunctionsTriangle;
    }
}

uint32_t ShapeFunctionIndexer::getNumOfInternalShapeFunctions(ElementType elementType) const
{
    if (elementType == ElementType_Parallelogram)
    {
        return m_numOfInternalShapeFunctionsQuadrilateral;
    }
    else
    {
        return m_numOfInternalShapeFunctionsTriangle;
    }
}

uint32_t ShapeFunctionIndexer::getNumOfShapeFunctions(ElementType elementType) const
{
    if (elementType == ElementType_Parallelogram)
    {
        return m_numOfShapeFunctionsQuadrilateral;
    }
    else
    {
        return m_numOfShapeFunctionsTriangle;
    }
}

ShapeFunctionDescriptor ShapeFunctionIndexer::getShapeFunctionDescriptor(ElementType elementType, uint32_t shapeFunctionIdx) const
{
    assert(shapeFunctionIdx < getNumOfShapeFunctions(elementType));
    if (shapeFunctionIdx < getNumOfNodalShapeFunctions(elementType))
    {
        return NodalShapeFunctionDescriptor(shapeFunctionIdx);
    }
    shapeFunctionIdx -= getNumOfNodalShapeFunctions(elementType);
    if (shapeFunctionIdx < getNumOfSideShapeFunctions(elementType))
    {
        return SideShapeFunctionDescriptor(shapeFunctionIdx / (m_p-1), shapeFunctionIdx % (m_p-1) + 2);
    }
    shapeFunctionIdx -= getNumOfSideShapeFunctions(elementType);
    if (elementType == ElementType_Parallelogram)
    {
        return getInternalShapeFunctionDescriptorQuadrilateral(shapeFunctionIdx);
    }
    else
    {
        return getInternalShapeFunctionDescriptorTriangle(shapeFunctionIdx);
    }
}

uint32_t ShapeFunctionIndexer::getShapeFunctionIndex(ElementType elementType, const ShapeFunctionDescriptor& descriptor) const
{
    if (elementType == ElementType_Parallelogram)
    {
        return getShapeFunctionIndexQuad(descriptor);
    }
    else
    {
        return getShapeFunctionIndexTri(descriptor);
    }
}

InternalShapeFunctionDescriptor ShapeFunctionIndexer::getInternalShapeFunctionDescriptor(ElementType elementType, uint32_t internalShapeFunctionIdx) const
{
    if (elementType == ElementType_Parallelogram)
    {
        return getInternalShapeFunctionDescriptorQuadrilateral(internalShapeFunctionIdx);
    }
    else
    {
        return getInternalShapeFunctionDescriptorTriangle(internalShapeFunctionIdx);
    }
}

uint32_t ShapeFunctionIndexer::getInternalShapeFunctionIndex(ElementType elementType, const InternalShapeFunctionDescriptor& descriptor) const
{
    if (elementType == ElementType_Parallelogram)
    {
        return getInternalShapeFunctionIndexQuadrilateral(descriptor);
    }
    else
    {
        return getInternalShapeFunctionIndexTriangle(descriptor);
    }
}

InternalShapeFunctionDescriptor ShapeFunctionIndexer::getInternalShapeFunctionDescriptorTriangle(uint32_t internalShapeFunctionIdx) const
{
    if (m_polynomialSpaceType == PolynomialSpaceType_Product)
    {
        return getInternalShapeFunctionDescriptorTriangleProduct(internalShapeFunctionIdx);
    }
    else
    {
        return getInternalShapeFunctionDescriptorTriangleTrunk(internalShapeFunctionIdx);
    }
}

InternalShapeFunctionDescriptor ShapeFunctionIndexer::getInternalShapeFunctionDescriptorTriangleTrunk(uint32_t internalShapeFunctionIdx) const
{
    uint32_t colStart = 0;
    for (int k = 0; k <= m_p - 3; k++)
    {
        const uint32_t colHeight = m_p-2-k;
        if (internalShapeFunctionIdx < colStart + colHeight)
        {
            return InternalShapeFunctionDescriptor(k, internalShapeFunctionIdx - colStart);
        }
        colStart += colHeight;
    }
    assert(false);
    return InternalShapeFunctionDescriptor{};
}

InternalShapeFunctionDescriptor ShapeFunctionIndexer::getInternalShapeFunctionDescriptorTriangleProduct(uint32_t internalShapeFunctionIdx) const
{
    return InternalShapeFunctionDescriptor(internalShapeFunctionIdx / (m_p-1), internalShapeFunctionIdx % (m_p-1));
}

InternalShapeFunctionDescriptor ShapeFunctionIndexer::getInternalShapeFunctionDescriptorQuadrilateral(uint32_t internalShapeFunctionIdx) const
{
    if (m_polynomialSpaceType == PolynomialSpaceType_Product)
    {
        return getInternalShapeFunctionDescriptorQuadrilateralProduct(internalShapeFunctionIdx);
    }
    else
    {
        return getInternalShapeFunctionDescriptorQuadrilateralTrunk(internalShapeFunctionIdx);
    }
}

InternalShapeFunctionDescriptor ShapeFunctionIndexer::getInternalShapeFunctionDescriptorQuadrilateralTrunk(uint32_t internalShapeFunctionIdx) const
{
    uint32_t colStart = 0;
    for (int k = 2; k <= m_p - 2; k++)
    {
        const uint32_t colHeight = m_p-1-k;
        if (internalShapeFunctionIdx < colStart + colHeight)
        {
            return InternalShapeFunctionDescriptor(k, internalShapeFunctionIdx - colStart + 2);
        }
        colStart += colHeight;
    }
    assert(false);
    return InternalShapeFunctionDescriptor{};
}

InternalShapeFunctionDescriptor ShapeFunctionIndexer::getInternalShapeFunctionDescriptorQuadrilateralProduct(uint32_t internalShapeFunctionIdx) const
{
    return InternalShapeFunctionDescriptor(internalShapeFunctionIdx / (m_p-1) + 2, internalShapeFunctionIdx % (m_p-1) + 2);
}

uint32_t ShapeFunctionIndexer::getShapeFunctionIndexQuad(const ShapeFunctionDescriptor& desc) const
{
    return std::visit([this](const auto& arg) { return this->getShapeFunctionIndexQuad(arg); }, desc);
}

uint32_t ShapeFunctionIndexer::getShapeFunctionIndexQuad(const NodalShapeFunctionDescriptor& desc) const
{
    return desc.nodeIdx;
}

uint32_t ShapeFunctionIndexer::getShapeFunctionIndexQuad(const SideShapeFunctionDescriptor& desc) const
{
    return 4 + desc.sideIdx * (m_p - 1) + desc.k - 2;
}

uint32_t ShapeFunctionIndexer::getShapeFunctionIndexQuad(const InternalShapeFunctionDescriptor& desc) const
{
    return 4 + 4 * (m_p - 1) + getInternalShapeFunctionIndexQuadrilateral(desc);
}

uint32_t ShapeFunctionIndexer::getShapeFunctionIndexTri(const ShapeFunctionDescriptor& desc) const
{
    return std::visit([this](const auto& arg) { return this->getShapeFunctionIndexTri(arg); }, desc);
}

uint32_t ShapeFunctionIndexer::getShapeFunctionIndexTri(const NodalShapeFunctionDescriptor& desc) const
{
    return desc.nodeIdx;
}

uint32_t ShapeFunctionIndexer::getShapeFunctionIndexTri(const SideShapeFunctionDescriptor& desc) const
{
    return 3 + desc.sideIdx * (m_p - 1) + desc.k - 2;
}

uint32_t ShapeFunctionIndexer::getShapeFunctionIndexTri(const InternalShapeFunctionDescriptor& desc) const
{
    return 3 + 3 * (m_p - 1) + getInternalShapeFunctionIndexTriangle(desc);
}

uint32_t ShapeFunctionIndexer::getInternalShapeFunctionIndexTriangle(const InternalShapeFunctionDescriptor& desc) const
{
    if (m_polynomialSpaceType == PolynomialSpaceType_Product)
    {
        return getInternalShapeFunctionIndexTriangleProduct(desc);
    }
    else
    {
        return getInternalShapeFunctionIndexTriangleTrunk(desc);
    }
}

uint32_t ShapeFunctionIndexer::getInternalShapeFunctionIndexTriangleTrunk(const InternalShapeFunctionDescriptor& desc) const
{
    assert(m_p >= 3);
    uint32_t res = 0;
    for (int i = 0; i < desc.k; i++)
    {
        res += m_p - 2 - i;
    }
    res += desc.l;
    return res;
}

uint32_t ShapeFunctionIndexer::getInternalShapeFunctionIndexTriangleProduct(const InternalShapeFunctionDescriptor& desc) const
{
    assert(m_p >= 2);
    return desc.k * (m_p-1) + desc.l;
}

uint32_t ShapeFunctionIndexer::getInternalShapeFunctionIndexQuadrilateral(const InternalShapeFunctionDescriptor& desc) const
{
    if (m_polynomialSpaceType == PolynomialSpaceType_Product)
    {
        return getInternalShapeFunctionIndexQuadrilateralProduct(desc);
    }
    else
    {
        return getInternalShapeFunctionIndexQuadrilateralTrunk(desc);
    }
}

uint32_t ShapeFunctionIndexer::getInternalShapeFunctionIndexQuadrilateralTrunk(const InternalShapeFunctionDescriptor& desc) const
{
    assert(m_p >= 4);
    uint32_t res = 0;
    for (int i = 0; i < desc.k - 2; i++)
    {
        res += m_p - 3 - i;
    }
    res += desc.l - 2;
    return res;
}

uint32_t ShapeFunctionIndexer::getInternalShapeFunctionIndexQuadrilateralProduct(const InternalShapeFunctionDescriptor& desc) const
{
    assert(m_p >= 2);
    return (desc.k - 2) * (m_p-1) + desc.l - 2;
}
} // namespace fem
