#include "fem/assembly/DiracLoadVector.hpp"

#include "fem/basis/BasisFunctionFactory.hpp"
#include "fem/basis/BasisFunctionIndexer.hpp"
#include "fem/basis/ShapeFunctionIndexer.hpp"

namespace fem
{
VectorXmpq assembleDiracLoadVector(const FemContext& ctx, const Vector2mpq& x_0, const ShapeFunctionFactory& shapeFunctionFactory)
{
    const BasisFunctionIndexer basisFunctionIndexer(ctx);
    const ShapeFunctionIndexer shapeFunctionIndexer(ctx.p, ctx.polynomialSpaceType);
    const uint32_t numOfBasisFunctions = basisFunctionIndexer.getNumOfBasisFunctions();
    VectorXmpq res(numOfBasisFunctions);
    const Mesh& mesh = *(ctx.mesh);
    const Mesh::ElementIndex elementIdx = getIndexOfElementContainingPoint(mesh, x_0);
    const Element& element = mesh.getElement(elementIdx);
    const ElementType elementType = element.getElementType();
    const AffineMap Finv = element.getReferenceElementMap().inverse();
    for (int shapeFnIdx = 0; shapeFnIdx < shapeFunctionIndexer.getNumOfShapeFunctions(elementType); shapeFnIdx++)
    {
        const auto desc = shapeFunctionIndexer.getShapeFunctionDescriptor(elementType, shapeFnIdx);
        const Polynomial2D& v = shapeFunctionFactory.getShapeFunction(elementType, desc);
        const uint32_t basisFnIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFnIdx);
        res(basisFnIdx) = v(Finv(x_0));
        if (const auto* descp = std::get_if<SideShapeFunctionDescriptor>(&desc))
        {
            const auto adjacentElementIdx = mesh.getIndexOfAdjacentElement(elementIdx, descp->sideIdx);
            if (adjacentElementIdx.has_value() && elementIdx < adjacentElementIdx.value() && descp->k % 2 != 0)
            {
                res(basisFnIdx) *= -1;
            }
        }
    }
    return res;
}

VectorXmpq assembleDiracLoadVector(const FemContext& ctx, const Vector2mpq& x_0)
{
    const BasisFunctionFactory basisFunctionFactory(ctx);
    const BasisFunctionIndexer basisFunctionIndexer(ctx);
    const uint32_t numOfBasisFunctions = basisFunctionIndexer.getNumOfBasisFunctions();
    VectorXmpq res(numOfBasisFunctions);
    const Mesh& mesh = *(ctx.mesh);
    const Mesh::ElementIndex elementIdx = getIndexOfElementContainingPoint(mesh, x_0);
    const Element& element = mesh.getElement(elementIdx);
    const AffineMap Finv = element.getReferenceElementMap().inverse();
    for (int shapeFunctionIdx = 0; shapeFunctionIdx < basisFunctionIndexer.getNumOfShapeFunctions(elementIdx); shapeFunctionIdx++)
    {
        const Polynomial2D v = basisFunctionFactory.getShapeFunction(elementIdx, shapeFunctionIdx);
        const uint32_t basisFunctionIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
        res(basisFunctionIdx) = v(Finv(x_0));
    }
    return res;
}
} // namespace fem
