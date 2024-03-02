#pragma once

#include <cstdint>

#include "fem/assembly/DiracLoadVector.hpp"
#include "fem/assembly/NeumannLoadVector.hpp"
#include "fem/basis/FemContext.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
VectorXmpq extractSubLoadVector(const FemContext& ctx, const VectorXmpq& loadVector, uint32_t p);
} // namespace fem
