#pragma once

#include <vector>

#include <gsl/gsl_math.h>

#include "fem/math/Function.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem::ut
{
inline const std::vector<BivariateFunction> referenceBoundaryFunctions = {
    /* g(x,y) = 1/(1+y^2) */
    [](const Vector2mpq& x) -> mpq_class
    {
        return 1 / (1 + x(1)*x(1));
    },

    /* g(x,y) = 1/sqrt(2) * (x-y)/(x^2+y^2) */
    [](const Vector2mpq& x) -> mpq_class
    {
        return (1/M_SQRT2) * (x(0) - x(1)) / (x(0)*x(0) + x(1)*x(1));
    }
};
} // namespace fem::ut
