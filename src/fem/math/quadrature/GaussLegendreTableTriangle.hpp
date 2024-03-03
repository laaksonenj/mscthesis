#pragma once

#include <cstdint>
#include <vector>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
class GaussLegendreTableTriangle
{
public:
    const auto& getWeights() const { return m_weights; }
    const auto& getAbscissas() const { return m_abscissas; }

protected:
    GaussLegendreTableTriangle() = default;
    ~GaussLegendreTableTriangle() = default;

protected:
    std::vector<mpq_class> m_weights;
    std::vector<Vector2mpq> m_abscissas;
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

inline const GaussLegendreTableTriangleQuadMapped glTableTriQuadMapped{100};
inline const GaussLegendreTableTriangleCrowdingFree glTableTriCrowdingFree{100};
inline const GaussLegendreTableTriangle& defaultGLTableTri = glTableTriCrowdingFree;
} // namespace fem
