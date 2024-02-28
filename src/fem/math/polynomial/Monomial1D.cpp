#include <sstream>

#include "fem/math/polynomial/Monomial1D.hpp"
#include "fem/multiprecision/Arithmetic.hpp"

namespace fem
{
mpq_class Monomial1D::operator()(const mpq_class& t) const
{
    return coefficient * pow(t, degree);
}

Monomial1D operator*(const Monomial1D& lhs, const Monomial1D& rhs)
{
    return Monomial1D(lhs.coefficient * rhs.coefficient, lhs.degree + rhs.degree);
}

Monomial1D operator-(const Monomial1D& op)
{
    return Monomial1D(-op.coefficient, op.degree);
}

bool operator==(const Monomial1D& lhs, const Monomial1D& rhs)
{
    const bool coefficientsAreEqual = lhs.coefficient == rhs.coefficient;
    const bool degreesAreEqual = lhs.degree == rhs.degree;
    const bool bothCoefficientsAreZero = lhs.coefficient == 0 && rhs.coefficient == 0;
    return (coefficientsAreEqual && degreesAreEqual) || bothCoefficientsAreZero;
}

bool operator!=(const Monomial1D& lhs, const Monomial1D& rhs)
{
    return !(lhs == rhs);
}

std::string toString(const Monomial1D& monomial)
{
    const mpq_class& coefficient = monomial.coefficient;
    if (coefficient == 0)
    {
        return "0";
    }
    std::stringstream ss;
    const uint32_t degree = monomial.degree;
    if (degree > 0)
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
    if (degree > 0)
    {
        ss << "t";
        if (degree > 1)
        {
            ss << "^" << degree;
        }
    }
    return ss.str();
}

std::ostream& operator<<(std::ostream& out, const Monomial1D& monomial)
{
    return out << toString(monomial);
}
} // namespace fem
