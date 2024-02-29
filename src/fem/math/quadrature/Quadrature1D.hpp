#pragma once

#include "fem/math/Function.hpp"
#include "fem/math/quadrature/GaussLegendreTable1D.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
mpq_class integrateGaussLegendre(const UnivariateFunction& f, const mpq_class& a, const mpq_class& b, const GaussLegendreTable1D& glTable = defaultGLTable1D);
} // namespace fem
