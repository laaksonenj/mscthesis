#pragma once

#include <ostream>

namespace fem
{
enum PolynomialSpaceType
{
    PolynomialSpaceType_Product,
    PolynomialSpaceType_Trunk
};

inline std::ostream& operator<<(std::ostream& out, PolynomialSpaceType polynomialSpaceType)
{
    if (polynomialSpaceType == PolynomialSpaceType_Product)
    {
        out << "product";
    }
    else
    {
        out << "trunk";
    }
    return out;
}
} // namespace fem
