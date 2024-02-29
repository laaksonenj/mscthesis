#pragma once

#include <cstdint>
#include <utility>

#include "fem/domain/Element.hpp"
#include "fem/math/Function.hpp"
#include "fem/math/quadrature/GaussLegendreTableQuadrilateral.hpp"
#include "fem/math/quadrature/GaussLegendreTableTriangle.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
mpq_class integrateGaussLegendre(const BivariateFunction& f, const Element& element);
mpq_class integrateGaussLegendre(const BivariateFunction& f, const Parallelogram& quad, const GaussLegendreTableQuadrilateral& glTable = defaultGLTableQuad);
mpq_class integrateGaussLegendre(const BivariateFunction& f, const Triangle& tri, const GaussLegendreTableTriangle& glTable = defaultGLTableTri);
} // namespace fem
