#include "fem/assembly/NeumannLoadVector.hpp"

#include <vector>

#include "fem/basis/BasisFunctionFactory.hpp"
#include "fem/basis/BasisFunctionIndexer.hpp"
#include "fem/basis/ShapeFunctionIndexer.hpp"
#include "fem/math/Quadrature.hpp"

namespace fem
{
VectorXmpq assembleNeumannLoadVector(const FemContext& ctx,
                                     const BivariateFunction& g,
                                     Mesh::ElementIndex elementIdx,
                                     Mesh::SideIndex localSideIdx,
                                     const ShapeFunctionFactory& shapeFunctionFactory)
{
    const Mesh& mesh = *ctx.mesh;
    assert(!mesh.getIndexOfAdjacentElement(elementIdx, localSideIdx).has_value());
    const BasisFunctionIndexer basisFunctionIndexer(ctx);
    const ShapeFunctionIndexer shapeFunctionIndexer(ctx.p, ctx.polynomialSpaceType);
    const uint32_t numOfBasisFunctions = basisFunctionIndexer.getNumOfBasisFunctions();
    VectorXmpq res(numOfBasisFunctions);

    const Element& element = mesh.getElement(elementIdx);
    const ElementType elementType = element.getElementType();
    const Side side = element.getSide(localSideIdx);
    const AffineMap Finv = element.getReferenceElementMap().inverse();
    auto r = [&side](const mpq_class& t) -> Vector2mpq
    {
        return (1-t)/2 * side.a + (1+t)/2 * side.b;
    };
    const Vector2mpq d = side.b - side.a;
    const mpq_class rGradNorm = mpq_class(sqrt(mpf_class(d(0)*d(0) + d(1)*d(1)))) / 2;

    std::vector<ShapeFunctionDescriptor> descs;
    descs.push_back(NodalShapeFunctionDescriptor(localSideIdx));
    descs.push_back(NodalShapeFunctionDescriptor((localSideIdx + 1) % element.getNumOfNodes()));
    for (int k = 2; k <= ctx.p; k++)
    {
        descs.push_back(SideShapeFunctionDescriptor(localSideIdx, k));
    }

    for (const auto& desc : descs)
    {
        const uint32_t shapeFnIdx = shapeFunctionIndexer.getShapeFunctionIndex(elementType, desc);
        const Polynomial2D& shapeFn = shapeFunctionFactory.getShapeFunction(elementType, desc);
        auto v = [&shapeFn, &Finv](const Vector2mpq& x) -> mpq_class
        {
            return shapeFn(Finv(x));
        };
        auto f = [&g, &v, &r](const mpq_class& t) -> mpq_class
        {
            return g(r(t)) * v(r(t));
        };
        const uint32_t basisFunctionIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFnIdx);
        res(basisFunctionIdx) = rGradNorm * integrateGaussLegendre(f, -1, 1);
    }

    return res;
}

VectorXmpq assembleNeumannLoadVector(const FemContext& ctx,
                                     const GradientFunction& grad_u,
                                     const ShapeFunctionFactory& shapeFunctionFactory)
{
    const BasisFunctionIndexer basisFunctionIndexer(ctx);
    VectorXmpq res(basisFunctionIndexer.getNumOfBasisFunctions());
    const Mesh& mesh = *(ctx.mesh);
    const auto meshBoundary = getMeshBoundary(mesh);
    #pragma omp parallel for
    for (int i = 0; i < meshBoundary.size(); i++)
    {
        const auto& [elementIdx, localSideIdx] = meshBoundary[i];
        const Element& element = mesh.getElement(elementIdx);
        const Side side = element.getSide(localSideIdx);
        const Vector2mpq n = calculateNormal(side);
        auto g = [&grad_u, &n](const Vector2mpq& x) -> mpq_class
        {
            return grad_u(x).dot(n);
        };
        const VectorXmpq res_i = assembleNeumannLoadVector(ctx, g, elementIdx, localSideIdx, shapeFunctionFactory);
        #pragma omp critical
        {
            res += res_i;
        }
    }
    return res;
}

VectorXmpq assembleNeumannLoadVector(const FemContext& ctx,
                                     const BivariateFunction& g,
                                     Mesh::ElementIndex elementIdx,
                                     Mesh::SideIndex localSideIdx)
{
    const Mesh& mesh = *ctx.mesh;
    assert(!mesh.getIndexOfAdjacentElement(elementIdx, localSideIdx).has_value());
    const BasisFunctionFactory basisFunctionFactory(ctx);
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
        const Polynomial2D shapeFn = basisFunctionFactory.getShapeFunction(elementIdx, shapeFunctionIdx);
        auto v = [&shapeFn, &Finv](const Vector2mpq& x) -> mpq_class
        {
            return shapeFn(Finv(x));
        };
        auto f = [&g, &v, &r](const mpq_class& t) -> mpq_class
        {
            return g(r(t)) * v(r(t));
        };
        const uint32_t basisFunctionIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
        res(basisFunctionIdx) = rGradNorm * integrateGaussLegendre(f, -1, 1);
    }
    return res;
}
} // namespace fem
