#pragma once

#include <string>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
std::string generateCommonHeader();
std::string namespaceBegin();
std::string namespaceEnd();
std::string generateVariableDefinition(const MatrixXmpq& matrix, const std::string& variableName);
std::string generateVariableDefinition(const VectorXmpq& vector, const std::string& variableName);
} // namespace fem
