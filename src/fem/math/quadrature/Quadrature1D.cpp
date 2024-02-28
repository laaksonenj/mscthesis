#include "fem/math/quadrature/Quadrature1D.hpp"

namespace fem
{
mpq_class integrateGaussLegendre(const UnivariateFunction& f, const mpq_class& a, const mpq_class& b, const GaussLegendreTable1D& glTable)
{
    mpq_class res = 0;
    auto g = [&f, &a, &b](const mpq_class& t) -> mpq_class
    {
        return f(((1-t)/2)*a + ((1+t)/2)*b);
    };
    for (const auto& [w, x] : glTable.getWeightAbscissaPairs())
    {
        res += w * g(x);
    }
    res *= (b-a)/2;
    return res;
}
} // namespace fem
