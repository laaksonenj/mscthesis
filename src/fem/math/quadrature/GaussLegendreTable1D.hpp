#pragma once

#include <cstdint>
#include <utility>
#include <vector>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
class GaussLegendreTable1D
{
public:
    explicit GaussLegendreTable1D(uint32_t n);

    const auto& getWeightAbscissaPairs() const { return m_weightAbscissaPairs; }

private:
    std::vector<std::pair<mpq_class, mpq_class>> m_weightAbscissaPairs;
};

inline GaussLegendreTable1D glTable1D100{100};
} // namespace fem
