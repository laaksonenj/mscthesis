#include <sstream>

#include "fem/math/polynomial/Monomial2D.hpp"
#include "fem/multiprecision/Arithmetic.hpp"

namespace fem
{
mpq_class Monomial2D::operator()(const Vector2mpq& p) const
{
    return coefficient * pow(p(0), degreeOfX) * pow(p(1), degreeOfY);
}

Monomial2D operator*(const Monomial2D& lhs, const Monomial2D& rhs)
{
    return Monomial2D(lhs.coefficient * rhs.coefficient, lhs.degreeOfX + rhs.degreeOfX, lhs.degreeOfY + rhs.degreeOfY);
}

Monomial2D operator-(const Monomial2D& op)
{
    return Monomial2D(-op.coefficient, op.degreeOfX, op.degreeOfY);
}

bool operator==(const Monomial2D& lhs, const Monomial2D& rhs)
{
    const bool coefficientsAreEqual = lhs.coefficient == rhs.coefficient;
    const bool degreesAreEqual = lhs.degreeOfX == rhs.degreeOfX && lhs.degreeOfY == rhs.degreeOfY;
    const bool bothCoefficientsAreZero = lhs.coefficient == 0 && rhs.coefficient == 0;
    return (coefficientsAreEqual && degreesAreEqual) || bothCoefficientsAreZero;
}

bool operator!=(const Monomial2D& lhs, const Monomial2D& rhs)
{
    return !(lhs == rhs);
}

std::string toString(const Monomial2D& monomial)
{
    const mpq_class& coefficient = monomial.coefficient;
    if (coefficient == 0)
    {
        return "0";
    }
    std::stringstream ss;
    const uint32_t degreeOfX = monomial.degreeOfX;
    const uint32_t degreeOfY = monomial.degreeOfY;
    if (degreeOfX > 0 || degreeOfY > 0)
    {
        if (coefficient == -1)
        {
            ss << "-";
        }
        else if (coefficient != 1)
        {
            ss << coefficient;
        }
    }
    else
    {
        ss << coefficient;
    }
    if (degreeOfX > 0)
    {
        ss << "x";
        if (degreeOfX > 1)
        {
            ss << "^" << degreeOfX;
        }
    }
    if (degreeOfY > 0)
    {
        ss << "y";
        if (degreeOfY > 1)
        {
            ss << "^" << degreeOfY;
        }
    }
    return ss.str();
}

std::ostream& operator<<(std::ostream& out, const Monomial2D& monomial)
{
    return out << toString(monomial);
}
} // namespace fem
