#pragma once

#include <cstdint>
#include <vector>

#include "fem/basis/BasisFunctionDescriptor.hpp"
#include "fem/basis/PolynomialSpaceType.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/basis/ShapeFunctionIndexer.hpp"
#include "fem/domain/Mesh.hpp"

namespace fem
{
class BasisFunctionIndexer
{
public:
    explicit BasisFunctionIndexer(const Mesh& mesh, uint32_t p, PolynomialSpaceType polynomialSpaceType);

    uint32_t getNumOfBasisFunctions() const { return m_numOfBasisFunctions; }
    uint32_t getNumOfShapeFunctions(Mesh::ElementIndex elementIdx) const;
    uint32_t getBasisFunctionIndex(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx) const;
    uint32_t getBasisFunctionIndex(const BasisFunctionDescriptor& descriptor) const;
    BasisFunctionDescriptor getBasisFunctionDescriptor(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx) const;
    BasisFunctionDescriptor getBasisFunctionDescriptor(uint32_t basisFunctionIndex) const;

    /* For testing & reference data generation purposes only */
    Polynomial2D getShapeFunction(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx) const;
    Polynomial2D getShapeFunctionDerivative(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx, char variable) const;

private:
    uint32_t getBasisFunctionIndexVisit(Mesh::ElementIndex elementIdx, const NodalShapeFunctionDescriptor& desc) const;
    uint32_t getBasisFunctionIndexVisit(Mesh::ElementIndex elementIdx, const SideShapeFunctionDescriptor& desc) const;
    uint32_t getBasisFunctionIndexVisit(Mesh::ElementIndex elementIdx, const InternalShapeFunctionDescriptor& desc) const;
    uint32_t getBasisFunctionIndexVisit(const NodalBasisFunctionDescriptor& desc) const;
    uint32_t getBasisFunctionIndexVisit(const SideBasisFunctionDescriptor& desc) const;
    uint32_t getBasisFunctionIndexVisit(const InternalBasisFunctionDescriptor& desc) const;

private:
    const Mesh& m_mesh;
    uint32_t m_p;
    PolynomialSpaceType m_polynomialSpaceType;

    uint32_t m_numOfNodalBasisFunctions;
    uint32_t m_numOfSideBasisFunctions;
    std::vector<uint32_t> m_accumulatedNumsOfInternalBasisFunctions;
    uint32_t m_numOfBasisFunctions;

    ShapeFunctionFactory m_shapeFunctionFactory;
    ShapeFunctionIndexer m_shapeFunctionIndexer;
};
} // namespace fem
