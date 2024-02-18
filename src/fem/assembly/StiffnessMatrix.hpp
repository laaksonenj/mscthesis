#pragma once

#include <cstdint>

#include "fem/basis/PolynomialSpaceType.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/domain/Mesh.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
MatrixXmpq assembleStiffnessMatrix(const Mesh& mesh, uint32_t p, PolynomialSpaceType polynomialSpaceType, const ShapeFunctionFactory& shapeFunctionFactory);
MatrixXmpq assembleStiffnessMatrix(const Mesh& mesh, uint32_t p, PolynomialSpaceType polynomialSpaceType);
} // namespace fem
