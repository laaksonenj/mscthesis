#pragma once

#include "fem/basis/FemContext.hpp"
#include "fem/basis/ShapeFunctionEvaluator.hpp"
#include "fem/domain/Mesh.hpp"
#include "fem/math/Function.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
mpq_class computeSquaredL2ErrorOverElement(const FemContext& ctx,
                                           const VectorXmpq& coeffs,
                                           const BivariateFunction& exact,
                                           Mesh::ElementIndex elementIdx,
                                           const Vector2mpq& x_0,
                                           const ShapeFunctionEvaluator& shapeFunctionEvaluator);
} // namespace fem
