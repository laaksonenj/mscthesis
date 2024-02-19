#pragma once

#include <vector>

#include <gsl/gsl_math.h>

#include "fem/math/Function.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem::ut
{
inline const std::vector<BivariateFunction> referenceBoundaryFunctions = {
    /* g(x,y) = 1/(1+y^2) */
    [](const Vector2mpq& p) -> mpq_class
    {
        return 1 / (1 + p(1)*p(1));
    },

    /* g(x,y) = 1/sqrt(2) * (x-y)/(x^2+y^2) */
    [](const Vector2mpq& p) -> mpq_class
    {
        return (1/M_SQRT2) * (p(0) - p(1)) / (p(0)*p(0) + p(1)*p(1));
    }
};
} // namespace fem::ut
