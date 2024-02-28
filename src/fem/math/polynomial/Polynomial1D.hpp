#pragma once

#include <cstdint>
#include <ostream>
#include <ranges>
#include <string>
#include <unordered_map>

#include "fem/math/polynomial/Monomial1D.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
class Polynomial1D
{
public:
    explicit Polynomial1D(const std::string& polynomialStr);
    Polynomial1D(const Monomial1D& monomial);
    Polynomial1D(const mpq_class& constant);
    Polynomial1D(int constant);
    Polynomial1D();

    auto getMonomials() const
    {
        return std::views::values(m_monomials);
    }

    mpq_class operator()(const mpq_class& t) const;

    Polynomial1D& operator+=(const Polynomial1D& rhs);
    Polynomial1D& operator-=(const Polynomial1D& rhs);
    Polynomial1D& operator*=(const Polynomial1D& rhs);

    friend Polynomial1D operator*(const Polynomial1D& lhs, const Polynomial1D& rhs);
    friend bool operator==(const Polynomial1D& lhs, const Polynomial1D& rhs);

private:
    void parsePolynomialString(const std::string& polynomialStr);
    void parseMonomialString(const std::string& monomialStr);

    void addMonomial(Monomial1D monomial);

    void doAddition(const Polynomial1D& rhs);
    void doSubtraction(const Polynomial1D& rhs);

private:
    std::unordered_map<uint32_t, Monomial1D> m_monomials;
};

Polynomial1D operator+(const Polynomial1D& lhs, const Polynomial1D& rhs);
Polynomial1D operator-(const Polynomial1D& lhs, const Polynomial1D& rhs);
Polynomial1D operator-(const Polynomial1D& op);
bool operator!=(const Polynomial1D& lhs, const Polynomial1D& rhs);
std::string toString(const Polynomial1D& polynomial);
std::ostream& operator<<(std::ostream& out, const Polynomial1D& polynomial);
} // namespace fem
