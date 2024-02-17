#include <gtest/gtest.h>

#include <cmath>

#include <gsl/gsl_math.h>

#include "fem/math/Quadrature.hpp"

namespace fem::ut
{
TEST(QuadratureTest, GaussLegendreUnivariate)
{
    auto f1 = [](const mpq_class& t) -> mpq_class
    {
        return (1-t)/(1+t*t);
    };
    EXPECT_DOUBLE_EQ(integrateGaussLegendre(f1, 0, 1, 14).get_d(), (M_PI - std::log(4))/4);

    auto f2 = [](const mpq_class& t) -> mpq_class
    {
        const mpq_class num = mpq_class("35/2")*t*t*t*t*t - mpq_class(25)*t*t*t + mpq_class("15/2")*t;
        const mpq_class den = 1 + t*t;
        return num/den;
    };
    EXPECT_NEAR(integrateGaussLegendre(f2, mpq_class("-2/7"), mpq_class("3/4"), 14).get_d(), 0.332078882023690280349671958279953, 1e-7);
}

TEST(QuadratureTest, GaussLegendreBivariateQuadrilateral)
{
    auto f1 = [](const Vector2mpq& p) -> mpq_class
    {
        const mpf_class r = sqrt(mpf_class(p(0)*p(0)+p(1)*p(1)));
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return mpq_class(std::log(r.get_d()));
        }
    };
    const Parallelogram quad1 = Parallelogram(Node(-1, -1), Node(1, -1), Node(1, 1), Node(-1, 1));
    EXPECT_NEAR(integrateGaussLegendre(f1, quad1, 42).get_d(), -1.47211, 0.01);

    auto f2 = [](const Vector2mpq& p) -> mpq_class
    {
        mpq_class res("4/7");
        for (int i = 0; i < 14; i++)
        {
            res *= p(0);
        }
        res *= p(1)*p(1)*p(1);
        res -= mpq_class("19/9")*p(0);
        return res;
    };
    const Parallelogram quad2 = Parallelogram(Node(1, 1), Node(2, 2), Node(2, 3), Node(1, 2));
    EXPECT_DOUBLE_EQ(integrateGaussLegendre(f2, quad2, 10).get_d(), mpq_class("41840537/2380").get_d());
}

TEST(QuadratureTest, GaussLegendreBivariateTriangle)
{
    auto f1 = [](const Vector2mpq& p) -> mpq_class
    {
        return p(0)*mpq_class(sqrt(mpf_class(1-p(1))));
    };
    const Triangle tri1 = Triangle(Vector2mpq(1, 0), Vector2mpq(0, 1), Vector2mpq(0, 0));
    EXPECT_NEAR(integrateGaussLegendre(f1, tri1, 10).get_d(), 1.0/7, 1e-7);

    auto f2 = [](const Vector2mpq& p) -> mpq_class
    {
        return mpq_class(sqrt(mpf_class(p(0)+p(1))));
    };
    EXPECT_NEAR(integrateGaussLegendre(f2, tri1, 14).get_d(), 2.0/5, 1e-7);

    auto f3 = [](const Vector2mpq& p) -> mpq_class
    {
        return mpq_class(std::exp(mpf_class(p(0)+p(1)).get_d()));
    };
    const Triangle tri3 = Triangle(Vector2mpq(0, 0), Vector2mpq(std::sqrt(3), 1), Vector2mpq(0, 2));
    EXPECT_NEAR(integrateGaussLegendre(f3, tri3, 10).get_d(), 9.763139379, 1e-7);

    auto f4 = [](const Vector2mpq& p) -> mpq_class
    {
        const mpq_class a = 2*p(0)+p(1);
        return a*a*a;
    };
    const Triangle tri4 = Triangle(Vector2mpq(-4, 1), Vector2mpq(-2.5, -3), Vector2mpq(-1, 1));
    EXPECT_NEAR(integrateGaussLegendre(f4, tri4, 10).get_d(), -1128, 1e-7);

    auto f5 = [](const Vector2mpq& p) -> mpq_class
    {
        const mpq_class a = p(1)+1;
        return p(0)*p(0)*a*a*a;
    };
    const Triangle tri5 = Triangle(Vector2mpq(-3, -2), Vector2mpq(5, -1), Vector2mpq(-2, 1));
    EXPECT_NEAR(integrateGaussLegendre(f5, tri5, 10).get_d(), 4301.0/420, 1e-7);
}
} // namespace fem::ut
