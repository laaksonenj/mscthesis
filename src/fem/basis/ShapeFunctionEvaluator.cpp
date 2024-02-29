#include "fem/basis/ShapeFunctionEvaluator.hpp"

#include <vector>

namespace fem
{
ShapeFunctionEvaluator::ShapeFunctionEvaluator(const ShapeFunctionFactory& shapeFunctionFactory)
    : m_shapeFunctionFactory(shapeFunctionFactory)
{
}

const mpq_class& ShapeFunctionEvaluator::evaluateShapeFunction(ElementType elementType, const ShapeFunctionDescriptor& descriptor, const Vector2mpq& x)
{
    auto& cache = m_cache[elementType];
    if (!cache.contains(descriptor))
    {
        cache.emplace(descriptor, PointEvalMap{});
    }
    if (!cache.at(descriptor).contains(x))
    {
        const Polynomial2D& shapeFn = m_shapeFunctionFactory.getShapeFunction(elementType, descriptor);
        cache.at(descriptor).emplace(x, shapeFn(x));
    }
    return cache.at(descriptor).at(x);
}

void ShapeFunctionEvaluator::preEvaluateShapeFunctions(ElementType elementType, const std::vector<Vector2mpq>& points)
{
    auto& cache = m_cache[elementType];
    std::vector<ShapeFunctionDescriptor> descs;
    for (const auto& [desc, shapeFn] : m_shapeFunctionFactory.getShapeFunctions(elementType))
    {
        cache.emplace(desc, PointEvalMap{});
        descs.push_back(desc);
    }
    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < descs.size(); i++)
    {
        const auto& desc = descs[i];
        const Polynomial2D& shapeFn = m_shapeFunctionFactory.getShapeFunction(elementType, desc);
        auto& pointEvalMap = cache.at(desc);
        for (const auto& x : points)
        {
            pointEvalMap.emplace(x, shapeFn(x));
        }
    }
}
} // namespace fem
