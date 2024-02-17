#pragma once

#include "fem/multiprecision/Types.hpp"

namespace fem
{
struct AffineMap
{
    Matrix2mpq A;
    Vector2mpq b;

    Vector2mpq operator()(const Vector2mpq& p) const;
    AffineMap inverse() const;
};

AffineMap compose(const AffineMap& F, const AffineMap& G);
bool operator==(const AffineMap& lhs, const AffineMap& rhs);
bool operator!=(const AffineMap& lhs, const AffineMap& rhs);
} // namespace fem
