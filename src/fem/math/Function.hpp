#pragma once

#include <functional>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
using UnivariateFunction = std::function<mpq_class(const mpq_class&)>;
using BivariateFunction = std::function<mpq_class(const Vector2mpq&)>;
using GradientFunction = std::function<Vector2mpq(const Vector2mpq&)>;
} // namespace fem
