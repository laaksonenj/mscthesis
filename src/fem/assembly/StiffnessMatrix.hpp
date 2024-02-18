#pragma once

#include "fem/basis/BasisFunctionIndexer.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
MatrixXmpq assembleStiffnessMatrix(const BasisFunctionIndexer& basisFunctionIndexer, const ShapeFunctionFactory& shapeFunctionFactory);
MatrixXmpq assembleStiffnessMatrix(const BasisFunctionIndexer& basisFunctionIndexer);
} // namespace fem
