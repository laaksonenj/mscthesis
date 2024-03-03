#include "fem/basis/ShapeFunctionEvaluator.hpp"

#include <vector>

namespace fem
{
ShapeFunctionEvaluator::ShapeFunctionEvaluator(const ShapeFunctionFactory& shapeFunctionFactory)
    : m_shapeFunctionFactory(shapeFunctionFactory)
{
}

mpq_class ShapeFunctionEvaluator::evaluate(ElementType elementType, const ShapeFunctionDescriptor& descriptor, const Vector2mpq& x) const
{
    const auto& cache = m_cache[elementType];
    if (cache.contains(descriptor) && cache.at(descriptor).contains(x))
    {
        return cache.at(descriptor).at(x);
    }
    else
    {
        const Polynomial2D& shapeFn = m_shapeFunctionFactory.getShapeFunction(elementType, descriptor);
        return shapeFn(x);
    }
}

void ShapeFunctionEvaluator::preEvaluate(ElementType elementType, const std::vector<Vector2mpq>& points)
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
