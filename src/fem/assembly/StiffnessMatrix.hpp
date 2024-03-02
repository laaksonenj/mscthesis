#pragma once

#include <cstdint>

#include "fem/basis/FemContext.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
MatrixXmpq assembleStiffnessMatrix(const FemContext& ctx, const ShapeFunctionFactory& shapeFunctionFactory);
MatrixXmpq extractSubStiffnessMatrix(const FemContext& ctx, const MatrixXmpq& stiffnessMatrix, uint32_t p);

/* Slow reference implementation */
MatrixXmpq assembleStiffnessMatrix(const FemContext& ctx);
} // namespace fem
