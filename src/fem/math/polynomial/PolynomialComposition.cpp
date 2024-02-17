#include "fem/math/polynomial/PolynomialComposition.hpp"

#include <sstream>

namespace fem
{
namespace
{
template<typename T>
T compose(const Monomial1D& f, const T& g)
{
    T res(1);
    for (int i = 0; i < f.degree; i++)
    {
        res *= g;
    }
    res *= f.coefficient;
    return res;
}

Polynomial2D compose(const Monomial2D& f, const Polynomial2D& g_x, const Polynomial2D& g_y)
{
    Polynomial2D res(1);
    for (int i = 0; i < f.degreeOfX; i++)
    {
        res *= g_x;
    }
    for (int i = 0; i < f.degreeOfY; i++)
    {
        res *= g_y;
    }
    res *= f.coefficient;
    return res;
}

Polynomial1D compose(const Monomial2D& f, const Polynomial1D& g_x, const Polynomial1D& g_y)
{
    Polynomial1D res(1);
    for (int i = 0; i < f.degreeOfX; i++)
    {
        res *= g_x;
    }
    for (int i = 0; i < f.degreeOfY; i++)
    {
        res *= g_y;
    }
    res *= f.coefficient;
    return res;
}
} // namespace

Polynomial1D compose(const Polynomial1D& f, const Polynomial1D& g)
{
    Polynomial1D res;
    for (const auto& monomial : f.getMonomials())
    {
        res += compose(monomial, g);
    }
    return res;
}

Polynomial2D compose(const Polynomial1D& f, const Polynomial2D& g)
{
    Polynomial2D res;
    for (const auto& monomial : f.getMonomials())
    {
        res += compose(monomial, g);
    }
    return res;
}

Polynomial2D compose(const Polynomial2D& f, const Polynomial2D& g_x, const Polynomial2D& g_y)
{
    Polynomial2D res;
    for (const auto& monomial : f.getMonomials())
    {
        res += compose(monomial, g_x, g_y);
    }
    return res;
}

Polynomial1D compose(const Polynomial2D& f, const Polynomial1D& g_x, const Polynomial1D& g_y)
{
    Polynomial1D res;
    for (const auto& monomial : f.getMonomials())
    {
        res += compose(monomial, g_x, g_y);
    }
    return res;
}

Polynomial2D compose(const Polynomial2D& f, const AffineMap& G)
{
    std::stringstream g_x_str;
    g_x_str << G.A(0,0) << "x";
    if (G.A(0,1) >= 0)
    {
        g_x_str << "+";
    }
    g_x_str << G.A(0,1) << "y";
    if (G.b(0) >= 0)
    {
        g_x_str << "+";
    }
    g_x_str << G.b(0);
    const Polynomial2D g_x(g_x_str.str());
    
    std::stringstream g_y_str;
    g_y_str << G.A(1,0) << "x";
    if (G.A(1,1) >= 0)
    {
        g_y_str << "+";
    }
    g_y_str << G.A(1,1) << "y";
    if (G.b(1) >= 0)
    {
        g_y_str << "+";
    }
    g_y_str << G.b(1);
    const Polynomial2D g_y(g_y_str.str());

    return compose(f, g_x, g_y);
}
} // namespace fem
