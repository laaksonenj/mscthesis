#pragma once

#include "fem/basis/FemContext.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
VectorXmpq assembleDiracLoadVector(const FemContext& ctx, const Vector2mpq& x_0, const ShapeFunctionFactory& shapeFunctionFactory);
VectorXmpq assembleDiracLoadVector(const FemContext& ctx, const Vector2mpq& x_0);
} // namespace fem
