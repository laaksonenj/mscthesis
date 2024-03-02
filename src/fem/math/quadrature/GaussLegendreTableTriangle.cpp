#include "fem/math/quadrature/GaussLegendreTableTriangle.hpp"

#include <ranges>

#include "fem/math/quadrature/GaussLegendreTable1D.hpp"

namespace fem
{
GaussLegendreTableTriangleQuadMapped::GaussLegendreTableTriangleQuadMapped(uint32_t n)
{
    const GaussLegendreTable1D t(n);
    for (auto [wi, xi] : std::views::zip(t.getWeights(), t.getAbscissas()))
    {
        const mpq_class u = (1 + xi) / 2;
        for (auto [wj, xj] : std::views::zip(t.getWeights(), t.getAbscissas()))
        {
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
        const mpq_class& wi = t1.getWeights().at(i);
        const mpq_class& xi = t1.getAbscissas().at(i);
        const mpq_class u = (1 + xi) / 2;
        for (auto [wj, xj] : std::views::zip(t2.getWeights(), t2.getAbscissas()))
        {
            const mpq_class v = (1 - xi) * (1 + xj) / 4;
            const mpq_class w = (1 - xi) * wi * wj / 8;
            m_weights.push_back(w);
            m_abscissas.push_back(Vector2mpq{u, v});
        }
    }
}
} // namespace fem
