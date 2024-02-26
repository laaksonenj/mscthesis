#include <cassert>
#include <cmath>
#include <numbers>

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

double log(const mpq_class& x)
{
    // x = a/b
    // log(x) = log(a) - log(b)
    assert(x > 0);
    return log(x.get_num()) - log(x.get_den());
}

double log(const mpz_class& x)
{
    // x = d * 2^exp
    // log(x) = log(d) + exp * log(2)
    assert(x > 0);
    signed long int exp;
    const double d = mpz_get_d_2exp(&exp, x.get_mpz_t());
    return std::log(d) + exp * std::numbers::ln2;
}
} // namespace fem
