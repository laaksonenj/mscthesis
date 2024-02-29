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
    for (const auto elem : table.getWeightAbscissaPairs())
    {
        const mpq_class& w = std::get<0>(elem);
        const Vector2mpq& x = std::get<1>(elem);
        res += w * g(x);
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
