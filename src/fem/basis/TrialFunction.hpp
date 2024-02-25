#pragma once

#include "fem/basis/BasisFunctionIndexer.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/math/Function.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
mpq_class evaluateTrialFunction(const VectorXmpq& coefficients,
                                const BasisFunctionIndexer& basisFunctionIndexer,
                                const Vector2mpq& p,
                                const ShapeFunctionFactory& shapeFunctionFactory);

mpq_class integrateTrialFunction(const VectorXmpq& coefficients,
                                 const BasisFunctionIndexer& basisFunctionIndexer,
                                 const ShapeFunctionFactory& shapeFunctionFactory);

void normalizeTrialFunction(VectorXmpq& coefficients,
                            const BasisFunctionIndexer& basisFunctionIndexer,
                            const ShapeFunctionFactory& shapeFunctionFactory);

mpq_class calculateSquaredL2Error(const VectorXmpq& coefficients,
                                  const BivariateFunction& u,
                                  const BasisFunctionIndexer& basisFunctionIndexer,
                                  const ShapeFunctionFactory& shapeFunctionFactory);

mpq_class calculateSquaredL2Error(const VectorXmpq& coefficients,
                                  const BivariateFunction& u,
                                  Mesh::ElementIndex elementIdx,
                                  const BasisFunctionIndexer& basisFunctionIndexer,
                                  const ShapeFunctionFactory& shapeFunctionFactory,
                                  uint32_t numOfGaussLegendrePoints = 20);

/* For testing purposes / sanity check */
mpq_class evaluateTrialFunction(const VectorXmpq& coefficients, const BasisFunctionIndexer& basisFunctionIndexer, const Vector2mpq& p);
mpq_class integrateTrialFunction(const VectorXmpq& coefficients, const BasisFunctionIndexer& basisFunctionIndexer);
} // namespace fem
