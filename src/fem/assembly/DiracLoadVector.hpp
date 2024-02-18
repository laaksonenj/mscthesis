#pragma once

#include <cstdint>

#include "fem/basis/PolynomialSpaceType.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/domain/Mesh.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
VectorXmpq assembleDiracLoadVector(const Mesh& mesh, uint32_t p, PolynomialSpaceType polynomialSpaceType, const Vector2mpq& x_0, const ShapeFunctionFactory& shapeFunctionFactory);
VectorXmpq assembleDiracLoadVector(const Mesh& mesh, uint32_t p, PolynomialSpaceType polynomialSpaceType, const Vector2mpq& x_0);
} // namespace fem
