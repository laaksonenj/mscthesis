#include "fem/math/quadrature/GaussLegendreTableQuadrilateral.hpp"

#include "fem/math/quadrature/GaussLegendreTable1D.hpp"

namespace fem
{
GaussLegendreTableQuadrilateral::GaussLegendreTableQuadrilateral(uint32_t n)
{
    const GaussLegendreTable1D t(n);
    for (int i = 0; i < n; i++)
    {
        const mpq_class& wi = t.getWeights().at(i);
        const mpq_class& xi = t.getAbscissas().at(i);
        for (int j = 0; j < n; j++)
        {
            const mpq_class& wj = t.getWeights().at(j);
            const mpq_class& xj = t.getAbscissas().at(j);
            m_weights.push_back(wi * wj);
            m_abscissas.push_back(Vector2mpq{xi, xj});
        }
    }
}
} // namespace fem
