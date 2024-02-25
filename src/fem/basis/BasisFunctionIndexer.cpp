#include "fem/basis/BasisFunctionIndexer.hpp"

#include <cassert>

namespace fem
{
auto calculateAccumulatedNumsOfInternalBasisFunctions(const Mesh& mesh, uint32_t p, PolynomialSpaceType polynomialSpaceType)
{
    std::vector<uint32_t> res{0};
    ShapeFunctionIndexer shapeFunctionIndexer(p, polynomialSpaceType);
    for (int elementIdx = 0; elementIdx < mesh.getNumOfElements(); elementIdx++)
    {
        const Element& element = mesh.getElement(elementIdx);
        res.push_back(res.back() + shapeFunctionIndexer.getNumOfInternalShapeFunctions(element.getElementType()));
    }
    return res;
}

BasisFunctionIndexer::BasisFunctionIndexer(const FemContext& ctx)
    : m_mesh(ctx.mesh)
    , m_p(ctx.p)
    , m_polynomialSpaceType(ctx.polynomialSpaceType)
    , m_numOfNodalBasisFunctions(m_mesh->getNumOfNodes())
    , m_numOfSideBasisFunctions(m_mesh->getNumOfSides() * (m_p-1))
    , m_accumulatedNumsOfInternalBasisFunctions(calculateAccumulatedNumsOfInternalBasisFunctions(*m_mesh, m_p, m_polynomialSpaceType))
    , m_numOfBasisFunctions(m_numOfNodalBasisFunctions + m_numOfSideBasisFunctions + m_accumulatedNumsOfInternalBasisFunctions.back())
    , m_shapeFunctionIndexer(ShapeFunctionIndexer(m_p, m_polynomialSpaceType))
{
}

uint32_t BasisFunctionIndexer::getNumOfElements() const
{
    return m_mesh->getNumOfElements();
}

uint32_t BasisFunctionIndexer::getNumOfShapeFunctions(Mesh::ElementIndex elementIdx) const
{
    assert(elementIdx < m_mesh->getNumOfElements());
    const Element& element = m_mesh->getElement(elementIdx);
    return m_shapeFunctionIndexer.getNumOfShapeFunctions(element.getElementType());
}

uint32_t BasisFunctionIndexer::getBasisFunctionIndex(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx) const
{
    assert(shapeFunctionIdx < getNumOfShapeFunctions(elementIdx));
    const Element& element = m_mesh->getElement(elementIdx);
    const ShapeFunctionDescriptor desc = m_shapeFunctionIndexer.getShapeFunctionDescriptor(element.getElementType(), shapeFunctionIdx);
    return std::visit([this, elementIdx](const auto& arg) { return this->getBasisFunctionIndexVisit(elementIdx, arg); }, desc);
}

uint32_t BasisFunctionIndexer::getBasisFunctionIndex(const BasisFunctionDescriptor& descriptor) const
{
    return std::visit([this](const auto& arg) { return this->getBasisFunctionIndexVisit(arg); }, descriptor);
}

BasisFunctionDescriptor BasisFunctionIndexer::getBasisFunctionDescriptor(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx) const
{
    return getBasisFunctionDescriptor(getBasisFunctionIndex(elementIdx, shapeFunctionIdx));
}

BasisFunctionDescriptor BasisFunctionIndexer::getBasisFunctionDescriptor(uint32_t basisFunctionIndex) const
{
    assert(basisFunctionIndex < getNumOfBasisFunctions());
    if (basisFunctionIndex < m_numOfNodalBasisFunctions)
    {
        return NodalBasisFunctionDescriptor(basisFunctionIndex);
    }
    basisFunctionIndex -= m_numOfNodalBasisFunctions;
    if (basisFunctionIndex < m_numOfSideBasisFunctions)
    {
        return SideBasisFunctionDescriptor(basisFunctionIndex / (m_p - 1), basisFunctionIndex % (m_p - 1) + 2);
    }
    basisFunctionIndex -= m_numOfSideBasisFunctions;
    for (int i = 1; i < m_accumulatedNumsOfInternalBasisFunctions.size(); i++)
    {
        if (basisFunctionIndex < m_accumulatedNumsOfInternalBasisFunctions.at(i))
        {
            const uint32_t elementIdx = i - 1;
            const uint32_t internalShapeFunctionIdx = basisFunctionIndex - m_accumulatedNumsOfInternalBasisFunctions.at(elementIdx);
            const Element& element = m_mesh->getElement(elementIdx);
            const InternalShapeFunctionDescriptor desc = m_shapeFunctionIndexer.getInternalShapeFunctionDescriptor(element.getElementType(), internalShapeFunctionIdx);
            return InternalBasisFunctionDescriptor(elementIdx, desc.k, desc.l);
        }
    }
    assert(false);
    return BasisFunctionDescriptor{};
}

Polynomial2D BasisFunctionIndexer::getShapeFunction(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx) const
{
    assert(shapeFunctionIdx < getNumOfShapeFunctions(elementIdx));
    const Element& element = m_mesh->getElement(elementIdx);
    const ShapeFunctionDescriptor desc = m_shapeFunctionIndexer.getShapeFunctionDescriptor(element.getElementType(), shapeFunctionIdx);
    const Polynomial2D& shapeFunction = m_shapeFunctionFactory.getShapeFunction(element.getElementType(), desc);
    if (const auto* descp = std::get_if<SideShapeFunctionDescriptor>(&desc))
    {
        const uint32_t localSideIdx = descp->sideIdx;
        const auto adjacentElementIdx = m_mesh->getIndexOfAdjacentElement(elementIdx, localSideIdx);
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

Polynomial2D BasisFunctionIndexer::getShapeFunctionDerivative(Mesh::ElementIndex elementIdx, uint32_t shapeFunctionIdx, char variable) const
{
    assert(shapeFunctionIdx < getNumOfShapeFunctions(elementIdx));
    assert(variable == 'x' || variable == 'y');
    return diff(getShapeFunction(elementIdx, shapeFunctionIdx), variable);
}

uint32_t BasisFunctionIndexer::getBasisFunctionIndexVisit(Mesh::ElementIndex elementIdx, const NodalShapeFunctionDescriptor& desc) const
{
    const uint32_t nodeIdx = m_mesh->getGlobalNodeIndex(elementIdx, desc.nodeIdx);
    return getBasisFunctionIndexVisit(NodalBasisFunctionDescriptor(nodeIdx));
}

uint32_t BasisFunctionIndexer::getBasisFunctionIndexVisit(Mesh::ElementIndex elementIdx, const SideShapeFunctionDescriptor& desc) const
{
    const uint32_t sideIdx = m_mesh->getGlobalSideIndex(elementIdx, desc.sideIdx);
    return getBasisFunctionIndexVisit(SideBasisFunctionDescriptor(sideIdx, desc.k));
}

uint32_t BasisFunctionIndexer::getBasisFunctionIndexVisit(Mesh::ElementIndex elementIdx, const InternalShapeFunctionDescriptor& desc) const
{
    return getBasisFunctionIndexVisit(InternalBasisFunctionDescriptor(elementIdx, desc.k, desc.l));
}

uint32_t BasisFunctionIndexer::getBasisFunctionIndexVisit(const NodalBasisFunctionDescriptor& desc) const
{
    return desc.nodeIdx;
}

uint32_t BasisFunctionIndexer::getBasisFunctionIndexVisit(const SideBasisFunctionDescriptor& desc) const
{
    uint32_t res = m_numOfNodalBasisFunctions;
    res += desc.sideIdx * (m_p-1);
    res += desc.k - 2;
    return res;
}

uint32_t BasisFunctionIndexer::getBasisFunctionIndexVisit(const InternalBasisFunctionDescriptor& desc) const
{
    uint32_t res = m_numOfNodalBasisFunctions + m_numOfSideBasisFunctions;
    res += m_accumulatedNumsOfInternalBasisFunctions.at(desc.elementIdx);
    const Element& element = m_mesh->getElement(desc.elementIdx);
    res += m_shapeFunctionIndexer.getInternalShapeFunctionIndex(element.getElementType(), InternalShapeFunctionDescriptor(desc.k, desc.l));
    return res;
}
} // namespace fem
