#pragma once

#include "fem/basis/FemContext.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/basis/ShapeFunctionIndexer.hpp"

namespace fem
{
/* Do not use this in "production". For validation purposes only. */
class BasisFunctionFactory
{
public:
    explicit BasisFunctionFactory(const FemContext& ctx);

    uint32_t getNumOfElements() const;
    uint32_t getNumOfShapeFunctions(Mesh::ElementIndex elementIdx) const;
    Polynomial2D getShapeFunction(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx) const;
    Polynomial2D getShapeFunctionDerivative(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx, char variable) const;

private:
    FemContext m_ctx;
    ShapeFunctionFactory m_shapeFunctionFactory;
    ShapeFunctionIndexer m_shapeFunctionIndexer;
};
} // namespace fem
