#pragma once

#include "fem/basis/FemContext.hpp"
#include "fem/basis/ShapeFunctionEvaluator.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
mpq_class evaluateTrialFunction(const VectorXmpq& coefficients, const FemContext& ctx, const Vector2mpq& x, ShapeFunctionEvaluator& shapeFunctionEvaluator);
mpq_class integrateTrialFunction(const VectorXmpq& coefficients, const FemContext& ctx, const ShapeFunctionFactory& shapeFunctionFactory);
void normalizeTrialFunction(VectorXmpq& coefficients, const FemContext& ctx, const ShapeFunctionFactory& shapeFunctionFactory);

/* For testing purposes / sanity check */
mpq_class evaluateTrialFunction(const VectorXmpq& coefficients, const FemContext& ctx, const Vector2mpq& x);
mpq_class integrateTrialFunction(const VectorXmpq& coefficients, const FemContext& ctx);
} // namespace fem
