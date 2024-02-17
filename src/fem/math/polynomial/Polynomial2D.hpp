#pragma once

#include <cstdint>
#include <ostream>
#include <ranges>
#include <string>
#include <unordered_map>
#include <utility>

#include <boost/functional/hash.hpp>

#include "fem/math/polynomial/Monomial2D.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
class Polynomial2D
{
public:
    explicit Polynomial2D(const std::string& polynomialStr);
    Polynomial2D(const Monomial2D& monomial);
    Polynomial2D(const mpq_class& constant);
    Polynomial2D(int constant);
    Polynomial2D();

    auto getMonomials() const
    {
        return std::views::values(m_monomials);
    }
    
    mpq_class operator()(const Vector2mpq& p) const;

    Polynomial2D& operator+=(const Polynomial2D& rhs);
    Polynomial2D& operator-=(const Polynomial2D& rhs);
    Polynomial2D& operator*=(const Polynomial2D& rhs);

    friend Polynomial2D operator*(const Polynomial2D& lhs, const Polynomial2D& rhs);
    friend bool operator==(const Polynomial2D& lhs, const Polynomial2D& rhs);

private:
    using DegreePair = std::pair<uint32_t, uint32_t>;
    using MonomialMap = std::unordered_map<DegreePair, Monomial2D, boost::hash<DegreePair>>;

private:
    void parsePolynomialString(const std::string& polynomialStr);
    void parseMonomialString(const std::string& monomialStr);

    void addMonomial(Monomial2D monomial);

    void doAddition(const Polynomial2D& rhs);
    void doSubtraction(const Polynomial2D& rhs);

private:
    MonomialMap m_monomials;
};

Polynomial2D operator+(const Polynomial2D& lhs, const Polynomial2D& rhs);
Polynomial2D operator-(const Polynomial2D& lhs, const Polynomial2D& rhs);
Polynomial2D operator-(const Polynomial2D& op);
bool operator!=(const Polynomial2D& lhs, const Polynomial2D& rhs);
std::string toString(const Polynomial2D& polynomial);
std::ostream& operator<<(std::ostream& out, const Polynomial2D& polynomial);
} // namespace fem