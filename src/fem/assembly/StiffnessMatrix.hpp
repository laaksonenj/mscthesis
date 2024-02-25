#pragma once

#include <cstdint>

#include "fem/basis/FemContext.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
MatrixXmpq assembleStiffnessMatrix(const FemContext& ctx, ShapeFunctionFactory& shapeFunctionFactory);
MatrixXmpq extractSubStiffnessMatrix(const MatrixXmpq& stiffnessMatrix, const FemContext& ctx, uint32_t p);

/* Slow reference implementation */
MatrixXmpq assembleStiffnessMatrix(const FemContext& ctx);
} // namespace fem
