#pragma once

#include <cassert>
#include <string>

#include "fem/domain/Node.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
struct Side
{
    Node a, b;

    Side(const Node& a, const Node& b)
        : a(a), b(b) { assert(a != b); }
};

Vector2mpq calculateNormal(const Side& side);
bool areParallel(const Side& side1, const Side& side2);
bool areIntersecting(const Side& side1, const Side& side2);
bool operator==(const Side& lhs, const Side& rhs);
bool operator!=(const Side& lhs, const Side& rhs);
} // namespace fem
