#include "fem/domain/Side.hpp"

#include <cassert>
#include <optional>
#include <sstream>
#include <utility>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
namespace
{
struct Line
{
    Vector2mpq orig;
    Vector2mpq dir;

    Line(const Vector2mpq& orig, const Vector2mpq& dir)
        : orig(orig), dir(dir) { assert(dir != Vector2mpq(0, 0)); }
};

bool areParallel(const Line& line1, const Line& line2)
{
    Vector3mpq line1Dir, line2Dir;
    line1Dir << line1.dir, 0;
    line2Dir << line2.dir, 0;
    return line1Dir.cross(line2Dir) == Vector3mpq(0, 0, 0);
}

std::optional<mpq_class> findIntersectionOfLineAndPoint(const Line& line, const Vector2mpq& point)
{
    assert(line.dir(0) != 0 || line.dir(1) != 0);
    const Vector2mpq rhs = point - line.orig;
    if (line.dir(0) == 0)
    {
        if (rhs(0) == 0)
        {
            const mpq_class t = rhs(1) / line.dir(1);
            return t;
        }
        else
        {
            return {};
        }
    }
    else if (line.dir(1) == 0)
    {
        if (rhs(1) == 0)
        {
            const mpq_class t = rhs(0) / line.dir(0);
            return t;
        }
        else
        {
            return {};
        }
    }
    else
    {
        const mpq_class t0 = rhs(0) / line.dir(0);
        const mpq_class t1 = rhs(1) / line.dir(1);
        if (t0 == t1)
        {
            return t0;
        }
        else
        {
            return {};
        }
    }
}
} // namespace

bool areParallel(const Side& side1, const Side& side2)
{
    const Vector2mpq& a1 = side1.a;
    const Vector2mpq& b1 = side1.b;
    const Vector2mpq& a2 = side2.a;
    const Vector2mpq& b2 = side2.b;
    return areParallel(Line(a1, b1 - a1), Line(a2, b2 - a2));
}

bool areIntersecting(const Side& side1, const Side& side2)
{
    const Vector2mpq& a1 = side1.a;
    const Vector2mpq& b1 = side1.b;
    const Vector2mpq& a2 = side2.a;
    const Vector2mpq& b2 = side2.b;
    if (!areParallel(side1, side2))
    {
        Matrix2mpq mat;
        mat.col(0) = b1 - a1;
        mat.col(1) = a2 - b2;
        assert(mat.determinant() != 0);
        const Vector2mpq x = mat.partialPivLu().solve(a2 - a1);
        return (x(0) >= 0 && x(0) <= 1 && x(1) >= 0 && x(1) <= 1);
    }
    else
    {
        const Line line1(a1, b1 - a1);
        const auto t = findIntersectionOfLineAndPoint(line1, a2);
        if (t)
        {
            const auto s = findIntersectionOfLineAndPoint(line1, b2);
            assert(s);
            return (*t >= 0 && *t <= 1) || (*s >= 0 && *s <= 1);
        }
        else
        {
            return false;
        }
    }
}

bool operator==(const Side& lhs, const Side& rhs)
{
    return (lhs.a == rhs.a && lhs.b == rhs.b) || (lhs.a == rhs.b && lhs.b == rhs.a);
}

bool operator!=(const Side& lhs, const Side& rhs)
{
    return !(lhs == rhs);
}
} // namespace fem
