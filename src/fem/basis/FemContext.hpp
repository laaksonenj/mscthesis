#pragma once

#include <cassert>
#include <cstdint>
#include <memory>

#include "fem/basis/PolynomialSpaceType.hpp"
#include "fem/domain/Mesh.hpp"

namespace fem
{
struct FemContext
{
    std::shared_ptr<Mesh> mesh;
    uint32_t p;
    PolynomialSpaceType polynomialSpaceType;

    FemContext(const std::shared_ptr<Mesh>& mesh, uint32_t p, PolynomialSpaceType polynomialSpaceType)
        : mesh(mesh)
        , p(p)
        , polynomialSpaceType(polynomialSpaceType)
    {
        assert(mesh != nullptr);
    }
};
} // namespace fem
