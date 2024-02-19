#pragma once

#include "fem/basis/BasisFunctionIndexer.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/math/Function.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
VectorXmpq assembleNeumannLoadVector(const BasisFunctionIndexer& basisFunctionIndexer,
                                     const BivariateFunction& g,
                                     Mesh::ElementIndex elementIdx,
                                     Mesh::SideIndex localSideIdx,
                                     const ShapeFunctionFactory& shapeFunctionFactory);
VectorXmpq assembleNeumannLoadVector(const BasisFunctionIndexer& basisFunctionIndexer,
                                     const BivariateFunction& g,
                                     Mesh::ElementIndex elementIdx,
                                     Mesh::SideIndex localSideIdx);
} // namespace fem
