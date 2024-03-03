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
    const uint32_t n = glTable.getAbscissas().size();
    for (int i = 0; i < n; i++)
    {
        const mpq_class& w = glTable.getWeights().at(i);
        const mpq_class& x = glTable.getAbscissas().at(i);
        res += w * g(x);
    }
    res *= (b-a)/2;
    return res;
}
} // namespace fem
