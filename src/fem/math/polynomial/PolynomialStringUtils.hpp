#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
std::vector<std::string> getMonomialStrings(std::string polynomialStr);
uint32_t getDegreeOfVariableFromMonomialString(const std::string& monomialStr, char variable);
mpq_class getCoefficientFromMonomialString(const std::string& monomialStr);
} // namespace fem
