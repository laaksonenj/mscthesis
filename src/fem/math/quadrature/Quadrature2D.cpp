#include "fem/math/quadrature/Quadrature2D.hpp"

namespace fem
{
namespace
{
template<typename T>
mpq_class doQuadrature(const BivariateFunction& f, const Element& element, const T& table)
{
    mpq_class res = 0;
    const AffineMap F = element.getReferenceElementMap();
    auto g = [&f, &F](const Vector2mpq& x) -> mpq_class
    {
        return f(F(x));
    };
    const auto& weights = table.getWeights();
    const auto& abscissas = table.getAbscissas();
    const uint32_t n = weights.size();
    #pragma omp parallel
    {
        mpq_class subRes = 0;
        #pragma omp for
        for (int i = 0; i < n; i++)
        {
            const auto& w = weights.at(i);
            const auto& x = abscissas.at(i);
            subRes += w * g(x);
        }
        #pragma omp critical
        {
            res += subRes;
        }
    }
    res *= F.A.determinant();
    return res;
}
} // namespace

mpq_class integrateGaussLegendre(const BivariateFunction& f, const Element& element)
{
    if (element.getElementType() == ElementType_Parallelogram)
    {
        return integrateGaussLegendre(f, static_cast<const Parallelogram&>(element));
    }
    else
    {
        return integrateGaussLegendre(f, static_cast<const Triangle&>(element));
    }
}

mpq_class integrateGaussLegendre(const BivariateFunction& f, const Parallelogram& quad, const GaussLegendreTableQuadrilateral& glTable)
{
    return doQuadrature(f, quad, glTable);
}

mpq_class integrateGaussLegendre(const BivariateFunction& f, const Triangle& tri, const GaussLegendreTableTriangle& glTable)
{
    return doQuadrature(f, tri, glTable);
}
} // namespace fem
