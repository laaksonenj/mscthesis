#pragma once

#include "fem/math/AffineMap.hpp"
#include "fem/math/polynomial/Polynomial1D.hpp"
#include "fem/math/polynomial/Polynomial2D.hpp"

namespace fem
{
Polynomial1D compose(const Polynomial1D& f, const Polynomial1D& g);
Polynomial2D compose(const Polynomial1D& f, const Polynomial2D& g);
Polynomial2D compose(const Polynomial2D& f, const Polynomial2D& g_x, const Polynomial2D& g_y);
Polynomial1D compose(const Polynomial2D& f, const Polynomial1D& g_x, const Polynomial1D& g_y);
Polynomial2D compose(const Polynomial2D& f, const AffineMap& G);
}
