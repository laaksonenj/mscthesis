#include "apps/dirac/utils/GreensFunction.hpp"

#include <numbers>

#include "fem/math/Quadrature.hpp"
#include "fem/multiprecision/Arithmetic.hpp"

namespace fem
{
BivariateFunction getGreensFunction(const Vector2mpq& x_0)
{
    return [x_0](const Vector2mpq& x) -> mpq_class
    {
        if (x == x_0)
        {
            return 0;
        }
        const Vector2mpq d = x - x_0;
        return (-std::numbers::inv_pi/2) * log(mpq_class(Vector2mpf(d(0), d(1)).norm()));
    };
}

BivariateFunction getNormalizedGreensFunction(const Vector2mpq& x_0, const Mesh& mesh)
{
    const auto G = getGreensFunction(x_0);
    mpq_class integralOfG = 0;
    for (int elementIdx = 0; elementIdx < mesh.getNumOfElements(); elementIdx++)
    {
        const Element& element = mesh.getElement(elementIdx);
        integralOfG += integrateGaussLegendre(G, element);
    }
    const mpq_class normalizationConst = integralOfG / calculateMeshArea(mesh);
    return [G, normalizationConst](const Vector2mpq& x) -> mpq_class
    {
        return G(x) - normalizationConst;
    };
}

GradientFunction getGreensFunctionGradient(const Vector2mpq& x_0)
{
    return [x_0](const Vector2mpq& x) -> Vector2mpq
    {
        if (x == x_0)
        {
            return {0,0};
        }
        const Vector2mpq d = x - x_0;
        return (-std::numbers::inv_pi/2) * d / d.squaredNorm();
    };
}
} // namespace fem
