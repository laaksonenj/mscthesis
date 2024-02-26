#include "fem/basis/BasisFunctionFactory.hpp"

#include <cassert>

namespace fem
{
BasisFunctionFactory::BasisFunctionFactory(const FemContext& ctx)
    : m_ctx(ctx)
    , m_shapeFunctionIndexer(ShapeFunctionIndexer(ctx.p, ctx.polynomialSpaceType))
{
    const Mesh& mesh = *m_ctx.mesh;
    if (mesh.containsQuadrilateral())
    {
        m_shapeFunctionFactory.createShapeFunctions(ElementType_Parallelogram, m_ctx.p);
    }
    if (mesh.containsTriangle())
    {
        m_shapeFunctionFactory.createShapeFunctions(ElementType_Triangle, m_ctx.p);
    }
}

uint32_t BasisFunctionFactory::getNumOfElements() const
{
    return m_ctx.mesh->getNumOfElements();
}

uint32_t BasisFunctionFactory::getNumOfShapeFunctions(Mesh::ElementIndex elementIdx) const
{
    assert(elementIdx < getNumOfElements());
    const Element& element = m_ctx.mesh->getElement(elementIdx);
    return m_shapeFunctionIndexer.getNumOfShapeFunctions(element.getElementType());
}

Polynomial2D BasisFunctionFactory::getShapeFunction(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx) const
{
    assert(shapeFunctionIdx < getNumOfShapeFunctions(elementIdx));
    const Mesh& mesh = *m_ctx.mesh;
    const Element& element = mesh.getElement(elementIdx);
    const ShapeFunctionDescriptor desc = m_shapeFunctionIndexer.getShapeFunctionDescriptor(element.getElementType(), shapeFunctionIdx);
    const Polynomial2D& shapeFunction = m_shapeFunctionFactory.getShapeFunction(element.getElementType(), desc);
    if (const auto* descp = std::get_if<SideShapeFunctionDescriptor>(&desc))
    {
        const uint32_t localSideIdx = descp->sideIdx;
        const auto adjacentElementIdx = mesh.getIndexOfAdjacentElement(elementIdx, localSideIdx);
        if (adjacentElementIdx.has_value() && elementIdx < adjacentElementIdx.value())
        {
            if (element.getElementType() == ElementType_Parallelogram)
            {
                if (localSideIdx == 0 || localSideIdx == 2)
                {
                    return compose(shapeFunction, Polynomial2D("-x"), Polynomial2D("y"));
                }
                else
                {
                    return compose(shapeFunction, Polynomial2D("x"), Polynomial2D("-y"));
                }
            }
            else
            {
                if (localSideIdx == 0)
                {
                    return compose(shapeFunction, Polynomial2D("-x"), Polynomial2D("2-y"));
                }
                else if (localSideIdx == 1)
                {
                    return compose(shapeFunction, Polynomial2D("-x"), Polynomial2D("-y"));
                }
                else
                {
                    return compose(shapeFunction, Polynomial2D("2-x"), Polynomial2D("-y"));
                }
            }
        }
        else
        {
            return shapeFunction;
        }
    }
    else
    {
        return shapeFunction;
    }
}

Polynomial2D BasisFunctionFactory::getShapeFunctionDerivative(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx, char variable) const
{
    assert(variable == 'x' || variable == 'y');
    return diff(getShapeFunction(elementIdx, shapeFunctionIdx), variable);
}
} // namespace fem
