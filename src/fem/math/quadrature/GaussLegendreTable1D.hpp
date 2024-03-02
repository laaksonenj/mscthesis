#pragma once

#include <cstdint>
#include <vector>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
class GaussLegendreTable1D
{
public:
    explicit GaussLegendreTable1D(uint32_t n);

    const auto& getWeights() const { return m_weights; }
    const auto& getAbscissas() const { return m_abscissas; }

private:
    std::vector<mpq_class> m_weights;
    std::vector<mpq_class> m_abscissas;
};

inline const GaussLegendreTable1D glTable1D100{100};
inline const GaussLegendreTable1D& defaultGLTable1D = glTable1D100;
} // namespace fem
