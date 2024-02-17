#include "fem/math/Quadrature.hpp"

#include <cassert>

#include <gsl/gsl_integration.h>

namespace fem
{
namespace
{
mpq_class integrateOverReferenceQuadrilateralGaussLegendre(const BivariateFunction& f, uint32_t n)
{
    mpq_class res = 0;
    gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(n);
    double xi, wi, xj, wj;
    for (int i = 0; i < n; i++)
    {
        gsl_integration_glfixed_point(-1, 1, i, &xi, &wi, t);
        for (int j = 0; j < n; j++)
        {
            gsl_integration_glfixed_point(-1, 1, j, &xj, &wj, t);
            res += (wi * wj) * f(Vector2mpq(xi, xj));
        }
    }
    gsl_integration_glfixed_table_free(t);
    return res;
}

mpq_class integrateOverReferenceTriangleGaussLegendre(const BivariateFunction& f, uint32_t n)
{
    auto g = [&f](const Vector2mpq& p) -> mpq_class
    {
        const mpq_class u = (1+p(0))/2;
        const mpq_class v = (1+p(1))/2;
        const mpq_class x = (1-v/2)*u;
        const mpq_class y = (1-u/2)*v;
        const mpq_class C = (1-(2+p(0)+p(1))/4)/4;
        return C * f(Vector2mpq(x, y));
    };
    return integrateOverReferenceQuadrilateralGaussLegendre(g, n);
}

mpq_class integrateOverParallelogramGaussLegendre(const BivariateFunction& f, const Parallelogram& quad, uint32_t n)
{
    const AffineMap F = quad.getReferenceElementMap();
    auto g = [&f, &F](const Vector2mpq& p) -> mpq_class
    {
        return f(F(p));
    };
    return F.A.determinant() * integrateOverReferenceQuadrilateralGaussLegendre(g, n);
}

mpq_class integrateOverTriangleGaussLegendre(const BivariateFunction& f, const Triangle& tri, uint32_t n)
{
    const AffineMap F = tri.getReferenceElementMap();
    auto g = [&f, &F](const Vector2mpq& p) -> mpq_class
    {
        return f(F(p));
    };
    return F.A.determinant() * integrateOverReferenceTriangleGaussLegendre(g, n);
}
} // namespace

mpq_class integrateGaussLegendre(const UnivariateFunction& f, const mpq_class& a, const mpq_class& b, uint32_t n)
{
    mpq_class res = 0;
    gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(n);
    double xi, wi;
    const double ad = a.get_d();
    const double bd = b.get_d();
    for (int i = 0; i < n; i++)
    {
        gsl_integration_glfixed_point(ad, bd, i, &xi, &wi, t);
        res += wi * f(xi);
    }
    gsl_integration_glfixed_table_free(t);
    return res;
}

mpq_class integrateGaussLegendre(const BivariateFunction& f, const Element& element, uint32_t n)
{
    if (element.getElementType() == ElementType_Parallelogram)
    {
        return integrateOverParallelogramGaussLegendre(f, static_cast<const Parallelogram&>(element), n);
    }
    else if (element.getElementType() == ElementType_Triangle)
    {
        return integrateOverTriangleGaussLegendre(f, static_cast<const Triangle&>(element), n);
    }
    else
    {
        assert(false);
        return 0;
    }
}
} // namespace fem
