#pragma once

#include <cstdint>
#include <vector>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
class GaussLegendreTableQuadrilateral
{
public:
    explicit GaussLegendreTableQuadrilateral(uint32_t n);

    const auto& getWeights() const { return m_weights; }
    const auto& getAbscissas() const { return m_abscissas; }

private:
    std::vector<mpq_class> m_weights;
    std::vector<Vector2mpq> m_abscissas;
};

inline const GaussLegendreTableQuadrilateral glTableQuad100{100};
inline const GaussLegendreTableQuadrilateral& defaultGLTableQuad = glTableQuad100;
} // namespace fem
