#pragma once

#include "fem/domain/Mesh.hpp"
#include "fem/math/Function.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
BivariateFunction getGreensFunction(const Vector2mpq& x_0);
BivariateFunction getNormalizedGreensFunction(const Vector2mpq& x_0, const Mesh& mesh);
GradientFunction getGreensFunctionGradient(const Vector2mpq& x_0);
} // namespace fem
