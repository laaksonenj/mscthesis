#include "fem/math/polynomial/Polynomial2D.hpp"

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
    assert(polynomialStr.find_first_not_of("0123456789+-xy^/ ") == std::string::npos);
}
} // namespace

Polynomial2D::Polynomial2D(const std::string& polynomialStr)
{
    parsePolynomialString(polynomialStr);
}

Polynomial2D::Polynomial2D(const Monomial2D& monomial)
{
    addMonomial(monomial);
}

Polynomial2D::Polynomial2D(const mpq_class& constant)
    : Polynomial2D(Monomial2D(constant, 0, 0))
{
}

Polynomial2D::Polynomial2D(int constant)
    : Polynomial2D(mpq_class(constant))
{
}

Polynomial2D::Polynomial2D()
    : Polynomial2D(0)
{
}

mpq_class Polynomial2D::operator()(const Vector2mpq& p) const
{
    mpq_class res = 0;
    for (const auto& monomial : getMonomials())
    {
        res += monomial(p);
    }
    return res;
}

Polynomial2D& Polynomial2D::operator+=(const Polynomial2D& rhs)
{
    if (this != &rhs)
    {
        doAddition(rhs);
    }
    else
    {
        doAddition(Polynomial2D(rhs));
    }
    return *this;
}

Polynomial2D& Polynomial2D::operator-=(const Polynomial2D& rhs)
{
    if (this != &rhs)
    {
        doSubtraction(rhs);
    }
    else
    {
        doSubtraction(Polynomial2D(rhs));
    }
    return *this;
}

Polynomial2D& Polynomial2D::operator*=(const Polynomial2D& rhs)
{
    *this = *this * rhs;
    return *this;
}

void Polynomial2D::parsePolynomialString(const std::string& polynomialStr)
{
    validatePolynomialString(polynomialStr);
    for (const auto& monomialStr : getMonomialStrings(polynomialStr))
    {
        parseMonomialString(monomialStr);
    }
}

void Polynomial2D::parseMonomialString(const std::string& monomialStr)
{
    const mpq_class coefficient = getCoefficientFromMonomialString(monomialStr);
    const uint32_t degreeOfX = getDegreeOfVariableFromMonomialString(monomialStr, 'x');
    const uint32_t degreeOfY = getDegreeOfVariableFromMonomialString(monomialStr, 'y');
    addMonomial(Monomial2D(coefficient, degreeOfX, degreeOfY));
}

void Polynomial2D::addMonomial(const Monomial2D& monomial)
{
    const DegreePair key = std::make_pair(monomial.degreeOfX, monomial.degreeOfY);
    if (!m_monomials.contains(key))
    {
        m_monomials.emplace(key, monomial);
    }
    else
    {
        m_monomials.at(key).coefficient += monomial.coefficient;
    }
    if (mpq_sgn(m_monomials.at(key).coefficient.get_mpq_t()) == 0)
    {
        m_monomials.erase(key);
    }
}

void Polynomial2D::doAddition(const Polynomial2D& rhs)
{
    for (const auto& monomial : rhs.getMonomials())
    {
        addMonomial(monomial);
    }
}

void Polynomial2D::doSubtraction(const Polynomial2D& rhs)
{
    for (const auto& monomial : rhs.getMonomials())
    {
        addMonomial(-monomial);
    }
}

Polynomial2D operator+(const Polynomial2D& lhs, const Polynomial2D& rhs)
{
    Polynomial2D res = lhs;
    res += rhs;
    return res;
}

Polynomial2D operator-(const Polynomial2D& lhs, const Polynomial2D& rhs)
{
    Polynomial2D res = lhs;
    res -= rhs;
    return res;
}

Polynomial2D operator*(const Polynomial2D& lhs, const Polynomial2D& rhs)
{
    Polynomial2D res;
    for (const auto& lhsMonomial : lhs.getMonomials())
    {
        for (const auto& rhsMonomial : rhs.getMonomials())
        {
            res.addMonomial(lhsMonomial * rhsMonomial);
        }
    }
    return res;
}

Polynomial2D operator-(const Polynomial2D& op)
{
    return Polynomial2D(0) - op;
}

bool operator==(const Polynomial2D& lhs, const Polynomial2D& rhs)
{
    return lhs.m_monomials == rhs.m_monomials;
}

bool operator!=(const Polynomial2D& lhs, const Polynomial2D& rhs)
{
    return !(lhs == rhs);
}

std::string toString(const Polynomial2D& polynomial)
{
    auto polynomialIt = polynomial.getMonomials();
    std::vector<Monomial2D> monomials(polynomialIt.begin(), polynomialIt.end());
    std::sort(monomials.begin(), monomials.end(),
        [](const Monomial2D& lhs, const Monomial2D& rhs)
        {
            const auto lhsPair = std::make_pair(lhs.degreeOfX, lhs.degreeOfY);
            const auto rhsPair = std::make_pair(rhs.degreeOfX, rhs.degreeOfY);
            return lhsPair > rhsPair;
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

std::ostream& operator<<(std::ostream& out, const Polynomial2D& polynomial)
{
    return out << toString(polynomial);
}
} // namespace fem
