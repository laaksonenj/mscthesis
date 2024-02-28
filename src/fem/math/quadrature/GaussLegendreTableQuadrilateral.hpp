#pragma once

#include <cstdint>
#include <utility>
#include <vector>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
class GaussLegendreTableQuadrilateral
{
public:
    GaussLegendreTableQuadrilateral(uint32_t n);

    const auto& getWeightAbscissaPairs() const { return m_weightAbscissaPairs; }

private:
    std::vector<std::pair<mpq_class, Vector2mpq>> m_weightAbscissaPairs;
};

inline GaussLegendreTableQuadrilateral glTableQuad100{100};
} // namespace fem
