#pragma once

#include <cstdint>
#include <memory>

#include "fem/basis/PolynomialSpaceType.hpp"
#include "fem/basis/ShapeFunctionDescriptor.hpp"
#include "fem/domain/Element.hpp"

namespace fem
{
class ShapeFunctionIndexer
{
public:
    explicit ShapeFunctionIndexer(uint32_t p, PolynomialSpaceType polynomialSpaceType);

    uint32_t getNumOfNodalShapeFunctions(ElementType elementType) const;
    uint32_t getNumOfSideShapeFunctions(ElementType elementType) const;
    uint32_t getNumOfInternalShapeFunctions(ElementType elementType) const;
    uint32_t getNumOfShapeFunctions(ElementType elementType) const;

    ShapeFunctionDescriptor getShapeFunctionDescriptor(ElementType elementType, uint32_t shapeFunctionIdx) const;
    InternalShapeFunctionDescriptor getInternalShapeFunctionDescriptor(ElementType elementType, uint32_t internalShapeFunctionIdx) const;
    uint32_t getInternalShapeFunctionIndex(ElementType elementType, const InternalShapeFunctionDescriptor& descriptor) const;

private:
    InternalShapeFunctionDescriptor getInternalShapeFunctionDescriptorTriangle(uint32_t internalShapeFunctionIdx) const;
    InternalShapeFunctionDescriptor getInternalShapeFunctionDescriptorTriangleTrunk(uint32_t internalShapeFunctionIdx) const;
    InternalShapeFunctionDescriptor getInternalShapeFunctionDescriptorTriangleProduct(uint32_t internalShapeFunctionIdx) const;

    InternalShapeFunctionDescriptor getInternalShapeFunctionDescriptorQuadrilateral(uint32_t internalShapeFunctionIdx) const;
    InternalShapeFunctionDescriptor getInternalShapeFunctionDescriptorQuadrilateralTrunk(uint32_t internalShapeFunctionIdx) const;
    InternalShapeFunctionDescriptor getInternalShapeFunctionDescriptorQuadrilateralProduct(uint32_t internalShapeFunctionIdx) const;

    uint32_t getInternalShapeFunctionIndexTriangle(const InternalShapeFunctionDescriptor& desc) const;
    uint32_t getInternalShapeFunctionIndexTriangleTrunk(const InternalShapeFunctionDescriptor& desc) const;
    uint32_t getInternalShapeFunctionIndexTriangleProduct(const InternalShapeFunctionDescriptor& desc) const;

    uint32_t getInternalShapeFunctionIndexQuadrilateral(const InternalShapeFunctionDescriptor& desc) const;
    uint32_t getInternalShapeFunctionIndexQuadrilateralTrunk(const InternalShapeFunctionDescriptor& desc) const;
    uint32_t getInternalShapeFunctionIndexQuadrilateralProduct(const InternalShapeFunctionDescriptor& desc) const;

private:
    uint32_t m_p;
    PolynomialSpaceType m_polynomialSpaceType;

    uint32_t m_numOfNodalShapeFunctionsQuadrilateral;
    uint32_t m_numOfSideShapeFunctionsQuadrilateral;
    uint32_t m_numOfInternalShapeFunctionsQuadrilateral;
    uint32_t m_numOfShapeFunctionsQuadrilateral;

    uint32_t m_numOfNodalShapeFunctionsTriangle;
    uint32_t m_numOfSideShapeFunctionsTriangle;
    uint32_t m_numOfInternalShapeFunctionsTriangle;
    uint32_t m_numOfShapeFunctionsTriangle;
};
} // namespace fem
