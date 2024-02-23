#include "fem/assembly/LoadVector.hpp"

namespace fem
{
VectorXmpq extractSubLoadVector(const VectorXmpq& loadVector, const BasisFunctionIndexer& superBasisFunctionIndexer, const BasisFunctionIndexer& subBasisFunctionIndexer)
{
    const uint32_t numOfBasisFunctions = subBasisFunctionIndexer.getNumOfBasisFunctions();
    VectorXmpq res(numOfBasisFunctions);
    for (int i = 0; i < numOfBasisFunctions; i++)
    {
        const BasisFunctionDescriptor desc = subBasisFunctionIndexer.getBasisFunctionDescriptor(i);
        const uint32_t ii = superBasisFunctionIndexer.getBasisFunctionIndex(desc);
        res(i) = loadVector(ii);
    }
    return res;
}
} // namespace fem
