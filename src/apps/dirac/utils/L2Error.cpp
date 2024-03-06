#include "apps/dirac/utils/L2Error.hpp"

#include "fem/basis/TrialFunction.hpp"
#include "fem/math/Quadrature.hpp"
#include "fem/multiprecision/Arithmetic.hpp"

namespace fem
{
namespace
{
const std::vector<Vector2mpq>& getAbscissas(ElementType elementType)
{
    if (elementType == ElementType_Parallelogram)
    {
        return defaultGLTableQuad.getAbscissas();
    }
    else
    {
        return defaultGLTableTri.getAbscissas();
    }
}
} // namespace

mpq_class computeSquaredL2ErrorOverElement(const FemContext& ctx,
                                           const VectorXmpq& coeffs,
                                           const BivariateFunction& exact,
                                           Mesh::ElementIndex elementIdx,
                                           const Vector2mpq& x_0,
                                           const ShapeFunctionEvaluator& shapeFunctionEvaluator)
{
    mpq_class res = 0;
    auto f = [&ctx, &coeffs, &exact, &shapeFunctionEvaluator](const Vector2mpq& x) -> mpq_class
    {
        const mpq_class approx = evaluateTrialFunction(ctx, coeffs, x, shapeFunctionEvaluator);
        return pow(exact(x) - approx, 2);
    };
    const Mesh& mesh = *ctx.mesh;
    const Element& element = mesh.getElement(elementIdx);
    const auto subElements = element.subdivide(x_0);
    for (const auto& elementPtr : subElements)
    {
        res += integrateGaussLegendre(f, *elementPtr);
    }
    return res;
}

std::vector<Vector2mpq> getShapeFunctionEvaluationPointsForL2Error(ElementType elementType, const Mesh& mesh, const Vector2mpq& x_0)
{
    const std::vector<Vector2mpq>& abscissas = getAbscissas(elementType);
    std::vector<Vector2mpq> points = abscissas;
    for (int elementIdx = 0; elementIdx < mesh.getNumOfElements(); elementIdx++)
    {
        const Element& element = mesh.getElement(elementIdx);
        if (element.getElementType() != elementType)
        {
            continue;
        }
        const auto subdivision = element.subdivide(x_0);
        if (subdivision.size() == 1)
        {
            continue;
        }
        const AffineMap Finv = element.getReferenceElementMap().inverse();
        for (const auto& subElement : subdivision)
        {
            const AffineMap G = subElement->getReferenceElementMap();
            for (const auto& x : abscissas)
            {
                points.push_back(Finv(G(x)));
            }
        }
    }
    return points;
}
} // namespace fem