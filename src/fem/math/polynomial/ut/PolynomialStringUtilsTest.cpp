#include <gtest/gtest.h>

#include "fem/math/polynomial/PolynomialStringUtils.hpp"

namespace fem::ut
{
TEST(PolynomialStringUtilsTest, GetMonomialStrings)
{
    EXPECT_EQ(getMonomialStrings("-1/4t^21+t^3-t+t+5/7-0"), std::vector<std::string>({"-1/4t^21", "+t^3", "-t", "+t", "+5/7", "-0"}));
    EXPECT_EQ(getMonomialStrings("t"), std::vector<std::string>({"t"}));
    EXPECT_EQ(getMonomialStrings("0"), std::vector<std::string>({"0"}));
}

TEST(PolynomialStringUtilsTest, GetDegreeOfVariableFromMonomialString)
{
    EXPECT_EQ(getDegreeOfVariableFromMonomialString("-1/4t^21", 't'), 21);
    EXPECT_EQ(getDegreeOfVariableFromMonomialString("12/7t^9", 't'), 9);
    EXPECT_EQ(getDegreeOfVariableFromMonomialString("-5/7x^2024y^9", 'x'), 2024);
    EXPECT_EQ(getDegreeOfVariableFromMonomialString("-5/7x^2024y^9", 'y'), 9);
    EXPECT_EQ(getDegreeOfVariableFromMonomialString("t", 't'), 1);
    EXPECT_EQ(getDegreeOfVariableFromMonomialString("1", 't'), 0);
    EXPECT_EQ(getDegreeOfVariableFromMonomialString("0", 't'), 0);
}

TEST(PolynomialStringUtilsTest, GetCoefficientFromMonomialString)
{
    EXPECT_EQ(getCoefficientFromMonomialString("-1/4t^21"), mpq_class("-1/4"));
    EXPECT_EQ(getCoefficientFromMonomialString("-5/7x^2024y^9"), mpq_class("-5/7"));
    EXPECT_EQ(getCoefficientFromMonomialString("-t^9"), mpq_class(-1));
    EXPECT_EQ(getCoefficientFromMonomialString("t^9"), mpq_class(1));
    EXPECT_EQ(getCoefficientFromMonomialString("4"), mpq_class(4));
    EXPECT_EQ(getCoefficientFromMonomialString("-4"), mpq_class(-4));
    EXPECT_EQ(getCoefficientFromMonomialString("0"), mpq_class(0));
}
} // namespace fem::ut
