#pragma once

#include <cstdint>
#include <utility>
#include <vector>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
class GaussLegendreTableTriangle
{
public:
    const auto& getWeightAbscissaPairs() const { return m_weightAbscissaPairs; }

protected:
    GaussLegendreTableTriangle() = default;
    ~GaussLegendreTableTriangle() = default;

protected:
    std::vector<std::pair<mpq_class, Vector2mpq>> m_weightAbscissaPairs;
};

class GaussLegendreTableTriangleQuadMapped : public GaussLegendreTableTriangle
{
public:
    explicit GaussLegendreTableTriangleQuadMapped(uint32_t n);
};

class GaussLegendreTableTriangleCrowdingFree : public GaussLegendreTableTriangle
{
public:
    explicit GaussLegendreTableTriangleCrowdingFree(uint32_t n);
};

inline GaussLegendreTableTriangleQuadMapped glTableTriQuadMapped100{100};
inline GaussLegendreTableTriangleCrowdingFree glTableTriCrowdingFree100{100};
} // namespace fem
