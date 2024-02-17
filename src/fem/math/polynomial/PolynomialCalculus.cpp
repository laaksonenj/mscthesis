#include "fem/math/polynomial/PolynomialCalculus.hpp"

#include <algorithm>
#include <cassert>

namespace fem
{
namespace
{
Monomial1D diff(const Monomial1D& monomial)
{
    const int degree = monomial.degree;
    return Monomial1D(degree * monomial.coefficient, std::max(0, degree - 1));
}

Monomial2D diff(const Monomial2D& monomial, char var)
{
    const int degreeOfX = monomial.degreeOfX;
    const int degreeOfY = monomial.degreeOfY;
    uint32_t newDegreeOfX = degreeOfX;
    uint32_t newDegreeOfY = degreeOfY;
    mpq_class newCoefficient = monomial.coefficient;
    if (var == 'x')
    {
        newCoefficient *= degreeOfX;
        newDegreeOfX = std::max(0, degreeOfX - 1);
    }
    else if (var == 'y')
    {
        newCoefficient *= degreeOfY;
        newDegreeOfY = std::max(0, degreeOfY - 1);
    }
    else
    {
        assert(false);
    }
    return Monomial2D(newCoefficient, newDegreeOfX, newDegreeOfY);
}

mpq_class integrateOverReferenceTriangle(const Monomial2D& monomial)
{
    const uint32_t degreeOfX = monomial.degreeOfX;
    const uint32_t degreeOfY = monomial.degreeOfY;
    const mpz_class num = mpz_class::factorial(degreeOfX) * mpz_class::factorial(degreeOfY);
    const mpz_class den = mpz_class::factorial(degreeOfX + degreeOfY + 2);
    mpq_class res(num, den);
    res.canonicalize();
    res *= monomial.coefficient;
    return res;
}

mpq_class integrateOverReferenceQuadrilateral(const Monomial2D& monomial)
{
    const uint32_t degreeOfX = monomial.degreeOfX;
    const uint32_t degreeOfY = monomial.degreeOfY;
    if (degreeOfX % 2 != 0 || degreeOfY % 2 != 0)
    {
        return 0;
    }
    else
    {
        mpq_class res(4, (degreeOfX + 1) * (degreeOfY + 1));
        res.canonicalize();
        res *= monomial.coefficient;
        return res;
    }
}

mpq_class integrateOverReferenceTriangle(const Polynomial2D& polynomial)
{
    mpq_class res = 0;
    for (const auto& monomial : polynomial.getMonomials())
    {
        res += integrateOverReferenceTriangle(monomial);
    }
    return res;
}

mpq_class integrateOverReferenceQuadrilateral(const Polynomial2D& polynomial)
{
    mpq_class res = 0;
    for (const auto& monomial : polynomial.getMonomials())
    {
        res += integrateOverReferenceQuadrilateral(monomial);
    }
    return res;
}
} // namespace

Polynomial1D diff(const Polynomial1D& polynomial)
{
    Polynomial1D res;
    for (const auto& monomial : polynomial.getMonomials())
    {
        res += diff(monomial);
    }
    return res;
}

Polynomial2D diff(const Polynomial2D& polynomial, char variable)
{
    Polynomial2D res;
    for (const auto& monomial : polynomial.getMonomials())
    {
        res += diff(monomial, variable);
    }
    return res;
}

mpq_class integrateOverReferenceElement(const Polynomial2D& polynomial, ElementType elementType)
{
    if (elementType == ElementType_Triangle)
    {
        return integrateOverReferenceTriangle(polynomial);
    }
    else if (elementType == ElementType_Parallelogram)
    {
        return integrateOverReferenceQuadrilateral(polynomial);
    }
    else
    {
        assert(false);
        return 0;
    }
}
} // namespace fem
