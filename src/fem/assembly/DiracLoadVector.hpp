#pragma once

#include "fem/basis/BasisFunctionIndexer.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
VectorXmpq assembleDiracLoadVector(const BasisFunctionIndexer& basisFunctionIndexer, const Vector2mpq& x_0, const ShapeFunctionFactory& shapeFunctionFactory);
VectorXmpq assembleDiracLoadVector(const BasisFunctionIndexer& basisFunctionIndexer, const Vector2mpq& x_0);
} // namespace fem
