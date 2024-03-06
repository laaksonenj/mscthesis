#pragma once

#include "fem/multiprecision/Types.hpp"

namespace fem
{
mpq_class pow(const mpq_class& base, int exp);
double log(const mpq_class& x);
double log(const mpz_class& x);
} // namespace fem
