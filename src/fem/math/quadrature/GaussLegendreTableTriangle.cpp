#include "fem/math/quadrature/GaussLegendreTableTriangle.hpp"

#include "fem/math/quadrature/GaussLegendreTable1D.hpp"

namespace fem
{
GaussLegendreTableTriangleQuadMapped::GaussLegendreTableTriangleQuadMapped(uint32_t n)
{
    const GaussLegendreTable1D t(n);
    for (const auto& [wi, xi] : t.getWeightAbscissaPairs())
    {
        const mpq_class u = (1 + xi) / 2;
        for (const auto& [wj, xj] : t.getWeightAbscissaPairs())
        {
            const mpq_class v = (1 - xi) * (1 + xj) / 4;
            const mpq_class w = (1 - xi) * wi * wj / 8;
            m_weightAbscissaPairs.emplace_back(w, Vector2mpq{u, v});
        }
    }
}

GaussLegendreTableTriangleCrowdingFree::GaussLegendreTableTriangleCrowdingFree(uint32_t n)
{
    const GaussLegendreTable1D t1(n);
    for (int i = 0; i < n; i++)
    {
        const GaussLegendreTable1D t2(n - i);
        const auto& [wi, xi] = t1.getWeightAbscissaPairs().at(i);
        const mpq_class u = (1 + xi) / 2;
        for (const auto& [wj, xj] : t2.getWeightAbscissaPairs())
        {
            const mpq_class v = (1 - xi) * (1 + xj) / 4;
            const mpq_class w = (1 - xi) * wi * wj / 8;
            m_weightAbscissaPairs.emplace_back(w, Vector2mpq{u, v});
        }
    }
}
} // namespace fem
