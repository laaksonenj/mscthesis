#pragma once

#include <cstdint>
#include <ostream>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
struct Monomial1D
{
    mpq_class coefficient;
    uint32_t degree;

    mpq_class operator()(const mpq_class& t) const;
};

Monomial1D operator*(const Monomial1D& lhs, const Monomial1D& rhs);
Monomial1D operator-(const Monomial1D& op);
bool operator==(const Monomial1D& lhs, const Monomial1D& rhs);
bool operator!=(const Monomial1D& lhs, const Monomial1D& rhs);
std::string toString(const Monomial1D& monomial);
std::ostream& operator<<(std::ostream& out, const Monomial1D& monomial);
} // namespace fem
