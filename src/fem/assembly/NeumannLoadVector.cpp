#include "fem/assembly/NeumannLoadVector.hpp"

#include "fem/basis/BasisFunctionIndexer.hpp"
#include "fem/math/Quadrature.hpp"

namespace fem
{
VectorXmpq assembleNeumannLoadVector(const FemContext& ctx,
                                     const BivariateFunction& g,
                                     Mesh::ElementIndex elementIdx,
                                     Mesh::SideIndex localSideIdx,
                                     const ShapeFunctionFactory& shapeFunctionFactory)
{
    return assembleNeumannLoadVector(ctx, g, elementIdx, localSideIdx);
}

VectorXmpq assembleNeumannLoadVector(const FemContext& ctx,
                                     const GradientFunction& grad_u,
                                     const ShapeFunctionFactory& shapeFunctionFactory)
{
    const BasisFunctionIndexer basisFunctionIndexer(ctx);
    VectorXmpq res(basisFunctionIndexer.getNumOfBasisFunctions());
    const Mesh& mesh = *(ctx.mesh);
    for (const auto& [elementIdx, localSideIdx] : getMeshBoundary(mesh))
    {
        const Element& element = mesh.getElement(elementIdx);
        const Side side = element.getSide(localSideIdx);
        const Vector2mpq n = calculateNormal(side);
        auto g = [&grad_u, &n](const Vector2mpq& x) -> mpq_class
        {
            return grad_u(x).dot(n);
        };
        res += assembleNeumannLoadVector(ctx, g, elementIdx, localSideIdx, shapeFunctionFactory);
    }
    return res;
}

VectorXmpq assembleNeumannLoadVector(const FemContext& ctx,
                                     const BivariateFunction& g,
                                     Mesh::ElementIndex elementIdx,
                                     Mesh::SideIndex localSideIdx)
{
    const Mesh& mesh = *(ctx.mesh);
    assert(!mesh.getIndexOfAdjacentElement(elementIdx, localSideIdx).has_value());
    const BasisFunctionIndexer basisFunctionIndexer(ctx);
    const uint32_t numOfBasisFunctions = basisFunctionIndexer.getNumOfBasisFunctions();
    VectorXmpq res(numOfBasisFunctions);
    const Element& element = mesh.getElement(elementIdx);
    const AffineMap Finv = element.getReferenceElementMap().inverse();
    const Side side = element.getSide(localSideIdx);
    auto r = [&side](const mpq_class& t) -> Vector2mpq
    {
        return (1-t)/2 * side.a + (1+t)/2 * side.b;
    };
    const Vector2mpq d = side.b - side.a;
    const mpq_class rGradNorm = mpq_class(sqrt(mpf_class(d(0)*d(0) + d(1)*d(1)))) / 2;
    for (int shapeFunctionIdx = 0; shapeFunctionIdx < basisFunctionIndexer.getNumOfShapeFunctions(elementIdx); shapeFunctionIdx++)
    {
        const Polynomial2D shapeFn = basisFunctionIndexer.getShapeFunction(elementIdx, shapeFunctionIdx);
        auto v = [&shapeFn, &Finv](const Vector2mpq& x) -> mpq_class
        {
            return shapeFn(Finv(x));
        };
        auto f = [&g, &v, &r](const mpq_class& t) -> mpq_class
        {
            return g(r(t)) * v(r(t));
        };
        const uint32_t basisFunctionIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
        res(basisFunctionIdx) = rGradNorm * integrateGaussLegendre(f, -1, 1, 100);
    }
    return res;
}
} // namespace fem
