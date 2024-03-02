#pragma once

#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>

#include "fem/basis/PolynomialSpaceType.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
class ShapeFunctionEvaluator
{
public:
    explicit ShapeFunctionEvaluator(const ShapeFunctionFactory& shapeFunctionFactory);
    ShapeFunctionEvaluator(ShapeFunctionFactory&& shapeFunctionFactory) = delete;

    mpq_class evaluate(ElementType elementType, const ShapeFunctionDescriptor& descriptor, const Vector2mpq& x) const;
    void preEvaluate(ElementType elementType, const std::vector<Vector2mpq>& points);

private:
    const ShapeFunctionFactory& m_shapeFunctionFactory;

    struct Vector2mpqCompare
    {
        bool operator()(const Vector2mpq& lhs, const Vector2mpq& rhs) const
        {
            return lhs(0) != rhs(0) ? lhs(0) < rhs(0) : lhs(1) < rhs(1);
        }
    };
    using PointEvalMap = std::map<Vector2mpq, mpq_class, Vector2mpqCompare>;
    std::unordered_map<ShapeFunctionDescriptor, PointEvalMap> m_cache[2];
};
} // namespace fem
