#pragma once

#include <cstdint>

#include "fem/domain/Element.hpp"
#include "fem/math/Function.hpp"

namespace fem
{
mpq_class integrateGaussLegendre(const UnivariateFunction& f, const mpq_class& a, const mpq_class& b, uint32_t n);
mpq_class integrateGaussLegendre(const BivariateFunction& f, const Element& element, uint32_t n);
} // namespace fem
