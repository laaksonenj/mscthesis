#include "apps/dirac/utils/L2Error.hpp"

#include "fem/basis/TrialFunction.hpp"
#include "fem/math/Quadrature.hpp"
#include "fem/multiprecision/Arithmetic.hpp"

namespace fem
{
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
} // namespace fem