#pragma once

#include <cstdint>
#include <ranges>
#include <vector>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
class GaussLegendreTableQuadrilateral
{
public:
    GaussLegendreTableQuadrilateral(uint32_t n);

    auto getWeightAbscissaPairs() const
    {
        return std::views::zip(m_weights, m_abscissas) | std::views::as_const;
    }

    const auto& getAbscissas() const { return m_abscissas; }

private:
    std::vector<mpq_class> m_weights;
    std::vector<Vector2mpq> m_abscissas;
};

inline const GaussLegendreTableQuadrilateral glTableQuad100{100};
inline const GaussLegendreTableQuadrilateral& defaultGLTableQuad = glTableQuad100;
} // namespace fem
