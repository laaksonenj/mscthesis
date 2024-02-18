#include "fem/assembly/DiracLoadVector.hpp"

namespace fem
{
VectorXmpq assembleDiracLoadVector(const BasisFunctionIndexer& basisFunctionIndexer, const Vector2mpq& x_0, const ShapeFunctionFactory& shapeFunctionFactory)
{
    return assembleDiracLoadVector(basisFunctionIndexer, x_0);
}

VectorXmpq assembleDiracLoadVector(const BasisFunctionIndexer& basisFunctionIndexer, const Vector2mpq& x_0)
{
    const uint32_t numOfBasisFunctions = basisFunctionIndexer.getNumOfBasisFunctions();
    VectorXmpq res(numOfBasisFunctions);
    const Mesh& mesh = basisFunctionIndexer.getMesh();
    const Mesh::ElementIndex elementIdx = getIndexOfElementContainingPoint(mesh, x_0);
    const Element& element = mesh.getElement(elementIdx);
    const AffineMap Finv = element.getReferenceElementMap().inverse();
    for (int shapeFunctionIdx = 0; shapeFunctionIdx < basisFunctionIndexer.getNumOfShapeFunctions(elementIdx); shapeFunctionIdx++)
    {
        const Polynomial2D f = basisFunctionIndexer.getShapeFunction(elementIdx, shapeFunctionIdx);
        const uint32_t basisFunctionIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
        res(basisFunctionIdx) = f(Finv(x_0));
    }
    return res;
}
} // namespace fem
