#pragma once

#include "fem/assembly/DiracLoadVector.hpp"
#include "fem/assembly/NeumannLoadVector.hpp"
#include "fem/basis/BasisFunctionIndexer.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
VectorXmpq extractSubLoadVector(const VectorXmpq& loadVector, const BasisFunctionIndexer& superBasisFunctionIndexer, const BasisFunctionIndexer& subBasisFunctionIndexer);
} // namespace fem
