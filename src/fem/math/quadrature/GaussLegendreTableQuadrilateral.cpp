#include "fem/math/quadrature/GaussLegendreTableQuadrilateral.hpp"

#include <ranges>

#include "fem/math/quadrature/GaussLegendreTable1D.hpp"

namespace fem
{
GaussLegendreTableQuadrilateral::GaussLegendreTableQuadrilateral(uint32_t n)
{
    const GaussLegendreTable1D t(n);
    for (auto [wi, xi] : std::views::zip(t.getWeights(), t.getAbscissas()))
    {
        for (auto [wj, xj] : std::views::zip(t.getWeights(), t.getAbscissas()))
        {
            m_weights.push_back(wi * wj);
            m_abscissas.push_back(Vector2mpq{xi, xj});
        }
    }
}
} // namespace fem
