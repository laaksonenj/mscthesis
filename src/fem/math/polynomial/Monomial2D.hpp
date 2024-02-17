#pragma once

#include <cstdint>
#include <ostream>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
struct Monomial2D
{
    mpq_class coefficient;
    uint32_t degreeOfX, degreeOfY;

    mpq_class operator()(const Vector2mpq& p) const;
};

Monomial2D operator*(const Monomial2D& lhs, const Monomial2D& rhs);
Monomial2D operator-(const Monomial2D& op);
bool operator==(const Monomial2D& lhs, const Monomial2D& rhs);
bool operator!=(const Monomial2D& lhs, const Monomial2D& rhs);
std::string toString(const Monomial2D& monomial);
std::ostream& operator<<(std::ostream& out, const Monomial2D& monomial);
} // namespace fem
