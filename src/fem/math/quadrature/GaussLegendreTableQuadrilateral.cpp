#include "fem/math/quadrature/GaussLegendreTableQuadrilateral.hpp"

#include "fem/math/quadrature/GaussLegendreTable1D.hpp"

namespace fem
{
GaussLegendreTableQuadrilateral::GaussLegendreTableQuadrilateral(uint32_t n)
{
    const GaussLegendreTable1D t(n);
    for (const auto elem1 : t.getWeightAbscissaPairs())
    {
        const mpq_class& wi = std::get<0>(elem1);
        const mpq_class& xi = std::get<1>(elem1);
        for (const auto elem2 : t.getWeightAbscissaPairs())
        {
            const mpq_class& wj = std::get<0>(elem2);
            const mpq_class& xj = std::get<1>(elem2);
            m_weights.push_back(wi * wj);
            m_abscissas.push_back(Vector2mpq{xi, xj});
        }
    }
}
} // namespace fem
