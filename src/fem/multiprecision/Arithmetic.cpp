#include <cassert>
#include <cmath>
#include <numbers>

#include "fem/multiprecision/Arithmetic.hpp"

namespace fem
{
mpq_class pow(const mpq_class& base, int exp)
{
    assert(!(base == 0 && exp < 0));
    mpq_class res;
    if (exp >= 0)
    {
        mpz_pow_ui(res.get_num_mpz_t(), base.get_num_mpz_t(), exp);
        mpz_pow_ui(res.get_den_mpz_t(), base.get_den_mpz_t(), exp);
    }
    else
    {
        exp *= -1;
        mpz_pow_ui(res.get_num_mpz_t(), base.get_den_mpz_t(), exp);
        mpz_pow_ui(res.get_den_mpz_t(), base.get_num_mpz_t(), exp);
    }
    res.canonicalize();
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
