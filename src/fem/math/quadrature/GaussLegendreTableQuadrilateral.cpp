#include "fem/math/quadrature/GaussLegendreTableQuadrilateral.hpp"

#include "fem/math/quadrature/GaussLegendreTable1D.hpp"

namespace fem
{
GaussLegendreTableQuadrilateral::GaussLegendreTableQuadrilateral(uint32_t n)
{
    const GaussLegendreTable1D t(n);
    for (const auto& [wi, xi] : t.getWeightAbscissaPairs())
    {
        for (const auto& [wj, xj] : t.getWeightAbscissaPairs())
        {
            m_weightAbscissaPairs.emplace_back(wi * wj, Vector2mpq{xi, xj});
        }
    }
}
} // namespace fem
