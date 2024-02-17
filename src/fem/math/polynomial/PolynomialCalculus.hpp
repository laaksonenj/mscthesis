#pragma once

#include "fem/domain/Element.hpp"
#include "fem/math/polynomial/Polynomial1D.hpp"
#include "fem/math/polynomial/Polynomial2D.hpp"

namespace fem
{
Polynomial1D diff(const Polynomial1D& polynomial);
Polynomial2D diff(const Polynomial2D& polynomial, char variable);
mpq_class integrateOverReferenceElement(const Polynomial2D& polynomial, ElementType elementType);
} // namespace fem
