#include "fem/math/polynomial/Polynomial1D.hpp"

#include <cassert>

#include <boost/algorithm/string/trim.hpp>

#include "fem/math/polynomial/PolynomialStringUtils.hpp"

namespace fem
{
namespace
{
void validatePolynomialString(const std::string& polynomialStr)
{
    assert(polynomialStr.length() > 0);
    assert(polynomialStr.find_first_not_of("0123456789+-t^/ ") == std::string::npos);
}
} // namespace

Polynomial1D::Polynomial1D(const std::string& polynomialStr)
{
    parsePolynomialString(polynomialStr);
}

Polynomial1D::Polynomial1D(const Monomial1D& monomial)
{
    addMonomial(monomial);
}

Polynomial1D::Polynomial1D(const mpq_class& constant)
    : Polynomial1D(Monomial1D(constant, 0))
{
}

Polynomial1D::Polynomial1D(int constant)
    : Polynomial1D(mpq_class(constant))
{
}

Polynomial1D::Polynomial1D()
    : Polynomial1D(0)
{
}

mpq_class Polynomial1D::operator()(const mpq_class& t) const
{
    mpq_class res = 0;
    for (const auto& monomial : getMonomials())
    {
        res += monomial(t);
    }
    return res;
}

Polynomial1D& Polynomial1D::operator+=(const Polynomial1D& rhs)
{
    if (this != &rhs)
    {
        doAddition(rhs);
    }
    else
    {
        doAddition(Polynomial1D(rhs));
    }
    return *this;
}

Polynomial1D& Polynomial1D::operator-=(const Polynomial1D& rhs)
{
    if (this != &rhs)
    {
        doSubtraction(rhs);
    }
    else
    {
        doSubtraction(Polynomial1D(rhs));
    }
    return *this;
}

Polynomial1D& Polynomial1D::operator*=(const Polynomial1D& rhs)
{
    *this = *this * rhs;
    return *this;
}

void Polynomial1D::parsePolynomialString(const std::string& polynomialStr)
{
    validatePolynomialString(polynomialStr);
    for (const auto& monomialStr : getMonomialStrings(polynomialStr))
    {
        parseMonomialString(monomialStr);
    }
}

void Polynomial1D::parseMonomialString(const std::string& monomialStr)
{
    const mpq_class coefficient = getCoefficientFromMonomialString(monomialStr);
    const uint32_t degree = getDegreeOfVariableFromMonomialString(monomialStr, 't');
    addMonomial(Monomial1D(coefficient, degree));
}

void Polynomial1D::addMonomial(Monomial1D monomial)
{
    auto node = m_monomials.extract(monomial.degree);
    if (node)
    {
        const mpq_class& oldCoefficient = node.mapped().coefficient;
        monomial.coefficient += oldCoefficient;
    }
    if (monomial.coefficient != 0)
    {
        m_monomials.emplace(monomial.degree, monomial);
    }
}

void Polynomial1D::doAddition(const Polynomial1D& rhs)
{
    for (const auto& monomial : rhs.getMonomials())
    {
        addMonomial(monomial);
    }
}

void Polynomial1D::doSubtraction(const Polynomial1D& rhs)
{
    for (const auto& monomial : rhs.getMonomials())
    {
        addMonomial(-monomial);
    }
}

Polynomial1D operator+(const Polynomial1D& lhs, const Polynomial1D& rhs)
{
    Polynomial1D res = lhs;
    res += rhs;
    return res;
}

Polynomial1D operator-(const Polynomial1D& lhs, const Polynomial1D& rhs)
{
    Polynomial1D res = lhs;
    res -= rhs;
    return res;
}

Polynomial1D operator*(const Polynomial1D& lhs, const Polynomial1D& rhs)
{
    Polynomial1D res;
    for (const auto& lhsMonomial : lhs.getMonomials())
    {
        for (const auto& rhsMonomial : rhs.getMonomials())
        {
            res.addMonomial(lhsMonomial * rhsMonomial);
        }
    }
    return res;
}

Polynomial1D operator-(const Polynomial1D& op)
{
    return Polynomial1D(0) - op;
}

bool operator==(const Polynomial1D& lhs, const Polynomial1D& rhs)
{
    return lhs.m_monomials == rhs.m_monomials;
}

bool operator!=(const Polynomial1D& lhs, const Polynomial1D& rhs)
{
    return !(lhs == rhs);
}

std::string toString(const Polynomial1D& polynomial)
{
    auto polynomialIt = polynomial.getMonomials();
    std::vector<Monomial1D> monomials(polynomialIt.begin(), polynomialIt.end());
    std::sort(monomials.begin(), monomials.end(),
        [](const Monomial1D& lhs, const Monomial1D& rhs)
        {
            return lhs.degree > rhs.degree;
        }
    );
    std::stringstream ss;
    for (const auto& monomial : monomials)
    {
        if (monomial.coefficient > 0)
        {
            ss << "+";
        }
        ss << monomial;
    }
    std::string res = ss.str();
    boost::trim_left_if(res, [](const char c) { return c == '+'; });
    if (res.empty())
    {
        res = "0";
    }
    return res;
}

std::ostream& operator<<(std::ostream& out, const Polynomial1D& polynomial)
{
    return out << toString(polynomial);
}
} // namespace fem
