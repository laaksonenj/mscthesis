#include <gtest/gtest.h>

#include <numbers>

#include "fem/math/polynomial/Polynomial2D.hpp"
#include "fem/math/quadrature/Quadrature2D.hpp"
#include "fem/multiprecision/Arithmetic.hpp"

namespace fem::ut
{
TEST(Quadrature2DTest, IntegrationOverParallelogram1)
{
    auto f1 = [](const Vector2mpq& x) -> mpq_class
    {
        if (x == Vector2mpq{0, 0})
        {
            return 0;
        }
        else
        {
            const mpf_class r = sqrt(mpf_class(x(0)*x(0)+x(1)*x(1)));
            return log(mpq_class(r));
        }
    };
    const Parallelogram quad1 = Parallelogram(Node(-1, -1), Node(1, -1), Node(1, 1), Node(-1, 1));
    EXPECT_NEAR(integrateGaussLegendre(f1, quad1).get_d(), -1.472112985, 1e-3);

    auto f2 = [](const Vector2mpq& x) -> mpq_class
    {
        mpq_class res("4/7");
        for (int i = 0; i < 14; i++)
        {
            res *= x(0);
        }
        res *= x(1)*x(1)*x(1);
        res -= mpq_class("19/9")*x(0);
        return res;
    };
    const Parallelogram quad2 = Parallelogram(Node(1, 1), Node(2, 2), Node(2, 3), Node(1, 2));
    EXPECT_DOUBLE_EQ(integrateGaussLegendre(f2, quad2).get_d(), mpq_class("41840537/2380").get_d());
}

TEST(Quadrature2DTest, IntegrationOverParallelogram2)
{
    auto f = [](const Vector2mpq& x) -> mpq_class
    {
        if (x == Vector2mpq{0, 0})
        {
            return 0;
        }
        else
        {
            const mpf_class r = sqrt(mpf_class(x(0)*x(0)+x(1)*x(1)));
            return log(mpq_class(r));
        }
    };
    const Parallelogram quad1 = Parallelogram(Node(-1, -1), Node(0, -1), Node(0, 0), Node(-1, 0));
    const Parallelogram quad2 = Parallelogram(Node(0, -1), Node(1, -1), Node(1, 0), Node(0, 0));
    const Parallelogram quad3 = Parallelogram(Node(-1, 0), Node(0, 0), Node(0, 1), Node(-1, 1));
    const Parallelogram quad4 = Parallelogram(Node(0, 0), Node(1, 0), Node(1, 1), Node(0, 1));
    const mpq_class res = integrateGaussLegendre(f, quad1)
                          + integrateGaussLegendre(f, quad2)
                          + integrateGaussLegendre(f, quad3)
                          + integrateGaussLegendre(f, quad4);
    EXPECT_NEAR(res.get_d(), -1.472112985, 1e-7); // compare to 1e-3 above in IntegrationOverParallelogram1 => more accurate by keeping the singularity at a node
}

TEST(Quadrature2DTest, IntegrationOverParallelogram3)
{
    const Parallelogram quad = Parallelogram(Node(-1, -1), Node(1, -1), Node(1, 1), Node(-1, 1));
    const auto v = mpq_class("21/104") * Polynomial2D("-75x + 455x^3 - 819x^5 + 429x^7") * Polynomial2D("1 - y^2");
    Vector2mpq x_0;
    const auto u = [&x_0](const Vector2mpq& x) -> mpq_class
    {
        if (x == x_0)
        {
            return 0;
        }
        else
        {
            const Vector2mpq d = x - x_0;
            const mpf_class r = sqrt(mpf_class(d(0)*d(0) + d(1)*d(1)));
            return (-std::numbers::inv_pi/2) * log(mpq_class(r));
        }
    };
    const auto f = [&u, &v](const Vector2mpq& x) -> mpq_class
    {
        return pow(u(x) - v(x), 2);
    };

    /* x_0 at node */
    x_0 = Vector2mpq{1, -1};
    EXPECT_NEAR(integrateGaussLegendre(f, quad).get_d(), 6.73582732080, 1e-8);

    /* x_0 on side */
    x_0 = Vector2mpq("-1/3", "1");
    const Parallelogram quad1s = Parallelogram(Node(-1, -1), Node("-1/3", "-1"), Node("-1/3", "1"), Node(-1, 1));
    const Parallelogram quad2s = Parallelogram(Node("-1/3", "-1"), Node(1, -1), Node(1, 1), Node("-1/3", "1"));
    const mpq_class sideRes = integrateGaussLegendre(f, quad1s) + integrateGaussLegendre(f, quad2s);
    EXPECT_NEAR(sideRes.get_d(), 6.2802483954283, 1e-8);

    /* x_0 inside */
    x_0 = Vector2mpq("-3/5", "-1/4");
    const Parallelogram quad1i = Parallelogram(Node(-1, -1), Node("-3/5", "-1"), Node("-3/5", "-1/4"), Node("-1", "-1/4"));
    const Parallelogram quad2i = Parallelogram(Node("-3/5", "-1/4"), Node("-3/5", "-1"), Node(1, -1), Node("1", "-1/4"));
    const Parallelogram quad3i = Parallelogram(Node("-3/5", "-1/4"), Node("1", "-1/4"), Node(1, 1), Node("-3/5", "1"));
    const Parallelogram quad4i = Parallelogram(Node(-1, 1), Node("-1", "-1/4"), Node("-3/5", "-1/4"), Node("-3/5", "1"));
    const mpq_class intRes = integrateGaussLegendre(f, quad1i)
                             + integrateGaussLegendre(f, quad2i)
                             + integrateGaussLegendre(f, quad3i)
                             + integrateGaussLegendre(f, quad4i);
    EXPECT_NEAR(intRes.get_d(), 5.9946842232642, 1e-8);
}

TEST(Quadrature2DTest, IntegrationOverTriangle1)
{
    auto f1 = [](const Vector2mpq& x) -> mpq_class
    {
        return x(0)*mpq_class(sqrt(mpf_class(1-x(1))));
    };
    const Triangle tri1 = Triangle(Vector2mpq(1, 0), Vector2mpq(0, 1), Vector2mpq(0, 0));
    EXPECT_NEAR(integrateGaussLegendre(f1, tri1).get_d(), 1.0/7, 1e-9);

    auto f2 = [](const Vector2mpq& x) -> mpq_class
    {
        return mpq_class(sqrt(mpf_class(x(0)+x(1))));
    };
    const Triangle tri2 = Triangle(Vector2mpq(1, 0), Vector2mpq(0, 1), Vector2mpq(0, 0));
    EXPECT_NEAR(integrateGaussLegendre(f2, tri2).get_d(), 2.0/5, 1e-9);

    auto f3 = [](const Vector2mpq& x) -> mpq_class
    {
        return mpq_class(std::exp(mpf_class(x(0)+x(1)).get_d()));
    };
    const Triangle tri3 = Triangle(Vector2mpq(0, 0), Vector2mpq(std::sqrt(3), 1), Vector2mpq(0, 2));
    EXPECT_NEAR(integrateGaussLegendre(f3, tri3).get_d(), 9.763139379, 1e-9);

    auto f4 = [](const Vector2mpq& x) -> mpq_class
    {
        const mpq_class a = 2*x(0)+x(1);
        return a*a*a;
    };
    const Triangle tri4 = Triangle(Vector2mpq(-4, 1), Vector2mpq(-2.5, -3), Vector2mpq(-1, 1));
    EXPECT_NEAR(integrateGaussLegendre(f4, tri4).get_d(), -1128, 1e-9);

    auto f5 = [](const Vector2mpq& x) -> mpq_class
    {
        const mpq_class a = x(1)+1;
        return x(0)*x(0)*a*a*a;
    };
    const Triangle tri5 = Triangle(Vector2mpq(-3, -2), Vector2mpq(5, -1), Vector2mpq(-2, 1));
    EXPECT_NEAR(integrateGaussLegendre(f5, tri5).get_d(), 4301.0/420, 1e-9);
}

TEST(Quadrature2DTest, IntegrationOverTriangle2)
{
    const Triangle tri = Triangle(Node(-1, -1), Node(1, -1), Node(-1, 1));
    const auto v = mpq_class("21/104") * Polynomial2D("-75x + 455x^3 - 819x^5 + 429x^7") * Polynomial2D("1 - y^2");
    Vector2mpq x_0;
    const auto u = [&x_0](const Vector2mpq& x) -> mpq_class
    {
        if (x == x_0)
        {
            return 0;
        }
        else
        {
            const Vector2mpq d = x - x_0;
            const mpf_class r = sqrt(mpf_class(d(0)*d(0) + d(1)*d(1)));
            return (-std::numbers::inv_pi/2) * log(mpq_class(r));
        }
    };
    const auto f = [&u, &v](const Vector2mpq& x) -> mpq_class
    {
        return pow(u(x) - v(x), 2);
    };

    /* x_0 at node 0 */
    x_0 = Vector2mpq{-1, -1};
    EXPECT_NEAR(integrateGaussLegendre(f, tri, glTableTriQuadMapped100).get_d(), 3.172693694155439, 1e-8);
    EXPECT_NEAR(integrateGaussLegendre(f, tri, glTableTriCrowdingFree100).get_d(), 3.172693694155439, 1e-8);

    /* x_0 at node 1 */
    x_0 = Vector2mpq{1, -1};
    EXPECT_NEAR(integrateGaussLegendre(f, tri, glTableTriQuadMapped100).get_d(), 3.499196781457323, 1e-8);
    EXPECT_NEAR(integrateGaussLegendre(f, tri, glTableTriCrowdingFree100).get_d(), 3.499196781457323, 1e-8);

    /* x_0 at node 2 */
    x_0 = Vector2mpq{-1, 1};
    EXPECT_NEAR(integrateGaussLegendre(f, tri, glTableTriQuadMapped100).get_d(), 3.20880247922330621291, 1e-8);
    EXPECT_NEAR(integrateGaussLegendre(f, tri, glTableTriCrowdingFree100).get_d(), 3.20880247922330621291, 1e-8);
}

TEST(Quadrature2DTest, IntegrationOverTriangle3)
{
    const Triangle tri = Triangle(Node(-1, -1), Node(1, -1), Node(-1, 1));
    const auto v = mpq_class("21/104") * Polynomial2D("-75x + 455x^3 - 819x^5 + 429x^7") * Polynomial2D("1 - y^2");
    Vector2mpq x_0;
    const auto u = [&x_0](const Vector2mpq& x) -> mpq_class
    {
        if (x == x_0)
        {
            return 0;
        }
        else
        {
            const Vector2mpq d = x - x_0;
            const mpf_class r = sqrt(mpf_class(d(0)*d(0) + d(1)*d(1)));
            return (-std::numbers::inv_pi/2) * log(mpq_class(r));
        }
    };
    const auto f = [&u, &v](const Vector2mpq& x) -> mpq_class
    {
        return pow(u(x) - v(x), 2);
    };

    /* x_0 on side */
    x_0 = Vector2mpq{0, 0};
    const Triangle tri1s = Triangle(Node(-1, -1), Node(1, -1), Node(0, 0));
    const Triangle tri2s = Triangle(Node(-1, 1), Node(-1, -1), Node(0, 0));
    const mpq_class resSide = integrateGaussLegendre(f, tri1s) + integrateGaussLegendre(f, tri2s);
    EXPECT_NEAR(resSide.get_d(), 3.073822175169817, 1e-8);

    /* x_0 inside */
    x_0 = Vector2mpq("-3/5", "-1/4");
    const Triangle tri1i = Triangle(Node(1, -1), Node("-3/5", "-1/4"), Node(-1, -1));
    const Triangle tri2i = Triangle(Node(-1, 1), Node("-3/5", "-1/4"), Node(1, -1));
    const Triangle tri3i = Triangle(Node(-1, -1), Node("-3/5", "-1/4"), Node(-1, 1));
    const mpq_class resInt = integrateGaussLegendre(f, tri1i) + integrateGaussLegendre(f, tri2i) + integrateGaussLegendre(f, tri3i);
    EXPECT_NEAR(resInt.get_d(), 2.9072571636138239, 1e-8);
}
} // namespace fem::ut
