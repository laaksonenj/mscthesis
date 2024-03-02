#pragma once

#include "fem/basis/FemContext.hpp"
#include "fem/basis/ShapeFunctionEvaluator.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
mpq_class evaluateTrialFunction(const FemContext& ctx, const VectorXmpq& coefficients, const Vector2mpq& x, const ShapeFunctionEvaluator& shapeFunctionEvaluator);
mpq_class integrateTrialFunction(const FemContext& ctx, const VectorXmpq& coefficients, const ShapeFunctionFactory& shapeFunctionFactory);
void normalizeTrialFunction(const FemContext& ctx, VectorXmpq& coefficients, const ShapeFunctionFactory& shapeFunctionFactory);

/* For testing purposes / sanity check */
mpq_class evaluateTrialFunction(const FemContext& ctx, const VectorXmpq& coefficients, const Vector2mpq& x);
mpq_class integrateTrialFunction(const FemContext& ctx, const VectorXmpq& coefficients);
} // namespace fem
