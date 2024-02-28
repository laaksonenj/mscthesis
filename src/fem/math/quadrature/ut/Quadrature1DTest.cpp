#include <gtest/gtest.h>

#include <numbers>

#include "fem/math/polynomial/Polynomial1D.hpp"
#include "fem/math/quadrature/Quadrature1D.hpp"
#include "fem/multiprecision/Arithmetic.hpp"

namespace fem::ut
{
TEST(Quadrature1DTest, GaussLegendre1)
{
    auto f1 = [](const mpq_class& t) -> mpq_class
    {
        return (1-t)/(1+t*t);
    };
    EXPECT_DOUBLE_EQ(integrateGaussLegendre(f1, 0, 1).get_d(), (std::numbers::pi - 2 * std::numbers::ln2)/4);

    auto f2 = [](const mpq_class& t) -> mpq_class
    {
        const mpq_class num = mpq_class("35/2")*t*t*t*t*t - mpq_class(25)*t*t*t + mpq_class("15/2")*t;
        const mpq_class den = 1 + t*t;
        return num/den;
    };
    EXPECT_NEAR(integrateGaussLegendre(f2, mpq_class("-2/7"), mpq_class("3/4")).get_d(), 0.332078882023690280349671958279953, 1e-10);
}

TEST(Quadrature1DTest, GaussLegendre2)
{
    const auto v = mpq_class("21/104") * Polynomial1D("-75t + 455t^3 - 819t^5 + 429t^7");
    const auto g = [](const mpq_class& t) -> mpq_class
    {
        return (-std::numbers::inv_pi/2) * 1 / (1 + pow((1-t)/8 + (1+t)/4, 2));
    };
    const auto f = [&v, &g](const mpq_class& t) -> mpq_class
    {
        return g(t) * v(t);
    };
    EXPECT_NEAR(integrateGaussLegendre(f, -1, 1).get_d(), -0.015248900493335452078, 1e-10);
}
} // namespace fem::ut
