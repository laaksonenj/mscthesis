#include "fem/assembly/LoadVector.hpp"

#include <cassert>

#include "fem/basis/BasisFunctionIndexer.hpp"

namespace fem
{
VectorXmpq extractSubLoadVector(const FemContext& ctx, const VectorXmpq& loadVector, uint32_t p)
{
    assert(p >= 1 && p <= ctx.p);
    const BasisFunctionIndexer superBasisFunctionIndexer(ctx);
    const BasisFunctionIndexer subBasisFunctionIndexer(FemContext(ctx.mesh, p, ctx.polynomialSpaceType));
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
