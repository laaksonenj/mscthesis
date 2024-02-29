#include "fem/math/quadrature/GaussLegendreTableTriangle.hpp"

#include "fem/math/quadrature/GaussLegendreTable1D.hpp"

namespace fem
{
GaussLegendreTableTriangleQuadMapped::GaussLegendreTableTriangleQuadMapped(uint32_t n)
{
    const GaussLegendreTable1D t(n);
    for (const auto elem1 : t.getWeightAbscissaPairs())
    {
        const mpq_class& wi = std::get<0>(elem1);
        const mpq_class& xi = std::get<1>(elem1);
        const mpq_class u = (1 + xi) / 2;
        for (const auto elem2 : t.getWeightAbscissaPairs())
        {
            const mpq_class& wj = std::get<0>(elem2);
            const mpq_class& xj = std::get<1>(elem2);
            const mpq_class v = (1 - xi) * (1 + xj) / 4;
            const mpq_class w = (1 - xi) * wi * wj / 8;
            m_weights.push_back(w);
            m_abscissas.push_back(Vector2mpq{u, v});
        }
    }
}

GaussLegendreTableTriangleCrowdingFree::GaussLegendreTableTriangleCrowdingFree(uint32_t n)
{
    const GaussLegendreTable1D t1(n);
    for (int i = 0; i < n; i++)
    {
        const GaussLegendreTable1D t2(n - i);
        const auto elem1 = t1.getWeightAbscissaPairs()[i];
        const mpq_class& wi = std::get<0>(elem1);
        const mpq_class& xi = std::get<1>(elem1);
        const mpq_class u = (1 + xi) / 2;
        for (const auto elem2 : t2.getWeightAbscissaPairs())
        {
            const mpq_class& wj = std::get<0>(elem2);
            const mpq_class& xj = std::get<1>(elem2);
            const mpq_class v = (1 - xi) * (1 + xj) / 4;
            const mpq_class w = (1 - xi) * wi * wj / 8;
            m_weights.push_back(w);
            m_abscissas.push_back(Vector2mpq{u, v});
        }
    }
}
} // namespace fem
