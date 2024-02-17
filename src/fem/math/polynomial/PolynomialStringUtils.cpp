#include "fem/math/polynomial/PolynomialStringUtils.hpp"

#include <cassert>
#include <algorithm>
#include <cctype>
#include <iostream>

#include <boost/algorithm/string/trim.hpp>

namespace fem
{
std::vector<std::string> getMonomialStrings(std::string polynomialStr)
{
    polynomialStr.erase(std::remove_if(polynomialStr.begin(),
                                       polynomialStr.end(),
                                       [](const char x) { return std::isspace(x); }),
                        polynomialStr.end());
    assert(polynomialStr.length() > 0);
    std::vector<std::string> monomialStrings;
    size_t pos;
    while ((pos = polynomialStr.find_first_of("+-", 1)) != std::string::npos)
    {
        std::string monomialStr = polynomialStr.substr(0, pos);
        monomialStrings.push_back(monomialStr);
        polynomialStr.erase(0, pos);
    }
    assert(polynomialStr.length() > 0);
    monomialStrings.push_back(polynomialStr);
    return monomialStrings;
}

uint32_t getDegreeOfVariableFromMonomialString(const std::string& monomialStr, char variable)
{
    const size_t varStartPos = monomialStr.find_first_of(variable);
    if (varStartPos == std::string::npos)
    {
        return 0;
    }
    const size_t expCharPos = varStartPos + 1;
    if (expCharPos == monomialStr.length() || monomialStr[expCharPos] != '^')
    {
        return 1;
    }
    const size_t degreeStartPos = expCharPos + 1;
    size_t degreeEndPos = monomialStr.find_first_not_of("0123456789", degreeStartPos);
    if (degreeEndPos == std::string::npos)
    {
        degreeEndPos = monomialStr.length();
    }
    const size_t degreeStrLen = degreeEndPos - degreeStartPos;
    const std::string degreeStr = monomialStr.substr(degreeStartPos, degreeStrLen);
    uint32_t degree;
    try
    {
        degree = std::stoul(degreeStr);
    }
    catch (...)
    {
        std::cout << "Invalid degree string " << degreeStr << std::endl;
        assert(false);
    }
    return degree;
}

mpq_class getCoefficientFromMonomialString(const std::string& monomialStr)
{
    assert(monomialStr.length() > 0);
    mpq_class res;
    std::string coeffStr;
    for (int i = 0; i < monomialStr.length(); i++)
    {
        if (!std::isalpha(monomialStr[i]))
        {
            coeffStr += monomialStr[i];
        }
        else
        {
            break;
        }
    }
    boost::trim_left_if(coeffStr, [](const char c) { return c == '+'; });
    if (coeffStr.empty())
    {
        res = mpq_class(1);
    }
    else if (coeffStr == "-")
    {
        res = mpq_class(-1);
    }
    else
    {
        try
        {
            res = mpq_class(coeffStr);
        }
        catch (...)
        {
            std::cout << "Invalid coefficient " << coeffStr << " in " << monomialStr << std::endl;
            assert(false);
        }
    }
    return res;
}
} // namespace fem
