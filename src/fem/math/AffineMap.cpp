#include "fem/math/AffineMap.hpp"

namespace fem
{
Vector2mpq AffineMap::operator()(const Vector2mpq& x) const
{
    const Vector2mpq Ax = A*x;
    return Ax + b;
}

AffineMap AffineMap::inverse() const
{
    const Matrix2mpq Ainv = A.inverse();
    return AffineMap(Ainv, -Ainv * b);
}

AffineMap compose(const AffineMap& F, const AffineMap& G)
{
    const Matrix2mpq ANew = F.A * G.A;
    Vector2mpq bNew = F.A * G.b;
    bNew += F.b;
    return AffineMap(ANew, bNew);
}

bool operator==(const AffineMap& lhs, const AffineMap& rhs)
{
    return (lhs.A == rhs.A) && (lhs.b == rhs.b);
}

bool operator!=(const AffineMap& lhs, const AffineMap& rhs)
{
    return !(lhs == rhs);
}
} // namespace fem
