#include "fem/basis/ShapeFunctionFactory.hpp"

#include <iostream>
#include <cassert>

namespace fem
{
ShapeFunctionFactory::ShapeFunctionFactory()
{
    m_legendrePolynomials.emplace(0, Polynomial1D(1));
    m_legendrePolynomials.emplace(1, Polynomial1D("t"));
}

const Polynomial2D& ShapeFunctionFactory::getShapeFunction(ElementType elementType, const ShapeFunctionDescriptor& descriptor) const
{
    if (elementType == ElementType_Parallelogram)
    {
        return getQuadShapeFunction(descriptor);
    }
    else
    {
        return getTriShapeFunction(descriptor);
    }
}

const Polynomial2D& ShapeFunctionFactory::getShapeFunctionDerivative(ElementType elementType, const ShapeFunctionDescriptor& descriptor, char variable) const
{
    assert(variable == 'x' || variable == 'y');
    if (elementType == ElementType_Parallelogram)
    {
        return getQuadShapeFunctionDerivative(descriptor, variable);
    }
    else
    {
        return getTriShapeFunctionDerivative(descriptor, variable);
    }
}

const Polynomial2D& ShapeFunctionFactory::getQuadShapeFunction(const ShapeFunctionDescriptor& d) const
{
    return std::visit([this](const auto& arg) -> const Polynomial2D& { return this->getQuadShapeFunction(arg); }, d);
}

const Polynomial2D& ShapeFunctionFactory::getQuadShapeFunction(const NodalShapeFunctionDescriptor& d) const
{
    const uint32_t nodeIdx = d.nodeIdx;
    assert(nodeIdx < 4);
    if (!m_quadShapeFunctions.contains(d))
    {
        if (nodeIdx == 0)
        {
            m_quadShapeFunctions.emplace(d, mpq_class("1/4") * Polynomial2D("1-x") * Polynomial2D("1-y"));
        }
        else if (nodeIdx == 1)
        {
            m_quadShapeFunctions.emplace(d, mpq_class("1/4") * Polynomial2D("1+x") * Polynomial2D("1-y"));
        }
        else if (nodeIdx == 2)
        {
            m_quadShapeFunctions.emplace(d, mpq_class("1/4") * Polynomial2D("1+x") * Polynomial2D("1+y"));
        }
        else
        {
            m_quadShapeFunctions.emplace(d, mpq_class("1/4") * Polynomial2D("1-x") * Polynomial2D("1+y"));
        }
    }
    return m_quadShapeFunctions.at(d);
}

const Polynomial2D& ShapeFunctionFactory::getQuadShapeFunction(const SideShapeFunctionDescriptor& d) const
{
    const uint32_t sideIdx = d.sideIdx;
    const uint32_t k = d.k;
    assert(sideIdx < 4);
    assert(k >= 2);
    if (!m_quadShapeFunctions.contains(d))
    {
        if (sideIdx == 0)
        {
            m_quadShapeFunctions.emplace(d, (mpq_class("1/2") * Polynomial2D("1-y")) * getPhi2D(k, 'x'));
        }
        else if (sideIdx == 1)
        {
            m_quadShapeFunctions.emplace(d, (mpq_class("1/2") * Polynomial2D("1+x")) * getPhi2D(k, 'y'));
        }
        else if (sideIdx == 2)
        {
            m_quadShapeFunctions.emplace(d, (mpq_class("1/2") * Polynomial2D("1+y")) * compose(getPhi(k), Polynomial2D("-x")));
        }
        else
        {
            m_quadShapeFunctions.emplace(d, (mpq_class("1/2") * Polynomial2D("1-x")) * compose(getPhi(k), Polynomial2D("-y")));
        }
    }
    return m_quadShapeFunctions.at(d);
}
    
const Polynomial2D& ShapeFunctionFactory::getQuadShapeFunction(const InternalShapeFunctionDescriptor& d) const
{
    const uint32_t k = d.k;
    const uint32_t l = d.l;
    assert(k >= 2);
    assert(l >= 2);
    if (!m_quadShapeFunctions.contains(d))
    {
        m_quadShapeFunctions.emplace(d, getPhi2D(k, 'x') * getPhi2D(l, 'y'));
    }
    return m_quadShapeFunctions.at(d);
}

const Polynomial2D& ShapeFunctionFactory::getTriShapeFunction(const ShapeFunctionDescriptor& d) const
{
    return std::visit([this](const auto& arg) -> const Polynomial2D& { return this->getTriShapeFunction(arg); }, d);
}

const Polynomial2D& ShapeFunctionFactory::getTriShapeFunction(const NodalShapeFunctionDescriptor& d) const
{
    const uint32_t nodeIdx = d.nodeIdx;
    assert(nodeIdx < 3);
    if (!m_triShapeFunctions.contains(d))
    {
        if (nodeIdx == 0)
        {
            m_triShapeFunctions.emplace(d, Polynomial2D("1-x-y"));
        }
        else if (nodeIdx == 1)
        {
            m_triShapeFunctions.emplace(d, Polynomial2D("x"));
        }
        else
        {
            m_triShapeFunctions.emplace(d, Polynomial2D("y"));
        }
    }
    return m_triShapeFunctions.at(d);
}

const Polynomial2D& ShapeFunctionFactory::getTriShapeFunction(const SideShapeFunctionDescriptor& d) const
{
    const uint32_t sideIdx = d.sideIdx;
    const uint32_t k = d.k;
    assert(sideIdx < 3);
    assert(k >= 2);
    if (!m_triShapeFunctions.contains(d))
    {
        const Polynomial1D& rho = getRho(k);
        const Polynomial2D& nodal1 = getTriShapeFunction(NodalShapeFunctionDescriptor(sideIdx));
        const Polynomial2D& nodal2 = getTriShapeFunction(NodalShapeFunctionDescriptor((sideIdx + 1) % 3));
        m_triShapeFunctions.emplace(d, (nodal1 * nodal2) * compose(rho, nodal2 - nodal1));
    }
    return m_triShapeFunctions.at(d);
}

const Polynomial2D& ShapeFunctionFactory::getTriShapeFunction(const InternalShapeFunctionDescriptor& d) const
{
    const uint32_t k = d.k;
    const uint32_t l = d.l;
    if (!m_triShapeFunctions.contains(d))
    {
        const Polynomial2D& nodal1 = getTriShapeFunction(NodalShapeFunctionDescriptor(0));
        const Polynomial2D& nodal2 = getTriShapeFunction(NodalShapeFunctionDescriptor(1));
        const Polynomial2D& nodal3 = getTriShapeFunction(NodalShapeFunctionDescriptor(2));
        const Polynomial2D& P_k = getShiftedLegendrePolynomial2D(k, 'x');
        const Polynomial2D& P_l = getShiftedLegendrePolynomial2D(l, 'y');
        m_triShapeFunctions.emplace(d, ((nodal1 * nodal2 * nodal3) * P_k) * P_l);
    }
    return m_triShapeFunctions.at(d);
}

const Polynomial2D& ShapeFunctionFactory::getQuadShapeFunctionDerivative(const ShapeFunctionDescriptor& d, char variable) const
{
    if (!m_quadShapeFunctionDerivatives.contains(d))
    {
        m_quadShapeFunctionDerivatives.emplace(d, std::unordered_map<char, Polynomial2D>{});
    }
    if (!m_quadShapeFunctionDerivatives.at(d).contains(variable))
    {
        m_quadShapeFunctionDerivatives.at(d).emplace(variable, diff(getQuadShapeFunction(d), variable));
    }
    return m_quadShapeFunctionDerivatives.at(d).at(variable);
}

const Polynomial2D& ShapeFunctionFactory::getTriShapeFunctionDerivative(const ShapeFunctionDescriptor& d, char variable) const
{
    if (!m_triShapeFunctionDerivatives.contains(d))
    {
        m_triShapeFunctionDerivatives.emplace(d, std::unordered_map<char, Polynomial2D>{});
    }
    if (!m_triShapeFunctionDerivatives.at(d).contains(variable))
    {
        m_triShapeFunctionDerivatives.at(d).emplace(variable, diff(getTriShapeFunction(d), variable));
    }
    return m_triShapeFunctionDerivatives.at(d).at(variable);
}

const Polynomial1D& ShapeFunctionFactory::getLegendrePolynomial(uint32_t n) const
{
    if (!m_legendrePolynomials.contains(n))
    {
        const mpq_class c1 = mpq_class(2*n - 1) / mpq_class(n);
        const mpq_class c2 = mpq_class(n - 1) / mpq_class(n);
        const Polynomial1D& PmMinus1 = getLegendrePolynomial(n-1);
        const Polynomial1D& PmMinus2 = getLegendrePolynomial(n-2);
        m_legendrePolynomials.emplace(n, (c1 * Polynomial1D("t")) * PmMinus1 - c2 * PmMinus2);
    }
    return m_legendrePolynomials.at(n);
}

const Polynomial2D& ShapeFunctionFactory::getShiftedLegendrePolynomial2D(uint32_t n, char variable) const
{
    assert(variable == 'x' || variable == 'y');
    if (!m_shiftedLegendrePolynomials2D.contains(n))
    {
        m_shiftedLegendrePolynomials2D.emplace(n, std::unordered_map<char, Polynomial2D>{});
    }
    if (!m_shiftedLegendrePolynomials2D.at(n).contains(variable))
    {
        const Polynomial2D g("2" + std::string{variable} + "-1");
        m_shiftedLegendrePolynomials2D.at(n).emplace(variable, compose(getLegendrePolynomial(n), g));
    }
    return m_shiftedLegendrePolynomials2D.at(n).at(variable);
}
    
const Polynomial1D& ShapeFunctionFactory::getPhi(uint32_t k) const
{
    assert(k >= 2);
    if (!m_phis.contains(k))
    {
        const mpq_class C = mpq_class((k-1) * k) / mpq_class(2*k - 1);
        m_phis.emplace(k, C * (getLegendrePolynomial(k) - getLegendrePolynomial(k-2)));
    }
    return m_phis.at(k);
}

const Polynomial2D& ShapeFunctionFactory::getPhi2D(uint32_t k, char variable) const
{
    assert(variable == 'x' || variable == 'y');
    assert(k >= 2);
    if (!m_phis2D.contains(k))
    {
        m_phis2D.emplace(k, std::unordered_map<char, Polynomial2D>{});
    }
    if (!m_phis2D.at(k).contains(variable))
    {
        m_phis2D.at(k).emplace(variable, compose(getPhi(k), Polynomial2D(std::string{variable})));
    }
    return m_phis2D.at(k).at(variable);
}

const Polynomial1D& ShapeFunctionFactory::getRho(uint32_t k) const
{
    assert(k >= 2);
    if (!m_rhos.contains(k))
    {
        m_rhos.emplace(k, -4 * diff(getLegendrePolynomial(k-1)));
    }
    return m_rhos.at(k);
}
} // namespace fem
