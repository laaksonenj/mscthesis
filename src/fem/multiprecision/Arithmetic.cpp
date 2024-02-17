#include <cassert>
#include <cmath>

#include "fem/multiprecision/Arithmetic.hpp"

namespace fem
{
mpq_class pow(mpq_class base, int exp)
{
    assert(!(base == 0 && exp < 0));
    mpq_class res = 1;
    if (exp < 0)
    {
        base = 1/base;
    }
    for (int i = 0; i < std::abs(exp); i++)
    {
        res *= base;
    }
    return res;
}
} // namespace fem
