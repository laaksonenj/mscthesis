#include "fem/basis/ShapeFunctionFactory.hpp"

#include <cassert>
#include <vector>

namespace fem
{
const Polynomial2D& ShapeFunctionFactory::getShapeFunction(ElementType elementType, const ShapeFunctionDescriptor& descriptor) const
{
    assert(m_shapeFunctions[elementType].contains(descriptor));
    return m_shapeFunctions[elementType].at(descriptor);
}

const Polynomial2D& ShapeFunctionFactory::getShapeFunctionDerivative(ElementType elementType, const ShapeFunctionDescriptor& descriptor, char variable) const
{
    assert(m_shapeFunctionDerivatives[elementType].contains(descriptor));
    assert(m_shapeFunctionDerivatives[elementType].at(descriptor).contains(variable));
    return m_shapeFunctionDerivatives[elementType].at(descriptor).at(variable);
}

void ShapeFunctionFactory::createShapeFunctions(ElementType elementType, uint32_t p)
{
    assert(p >= 1);
    if (elementType == ElementType_Parallelogram)
    {
        createQuadShapeFunctions(p);
    }
    else
    {
        createTriShapeFunctions(p);
    }

    /* Free up temporary memory */
    m_legendrePolynomials.clear();
    m_shiftedLegendrePolynomials2D.clear();
    m_phis.clear();
    m_phis2D.clear();
    m_rhos.clear();
}

void ShapeFunctionFactory::createQuadShapeFunctions(int p)
{
    createLegendrePolynomials(p);
    createPhis(p);

    std::vector<ShapeFunctionDescriptor> descs;
    for (int nodeIdx = 0; nodeIdx < 4; nodeIdx++)
    {
        descs.push_back(NodalShapeFunctionDescriptor(nodeIdx));
    }
    for (int sideIdx = 0; sideIdx < 4; sideIdx++)
    {
        for (int k = 2; k <= p; k++)
        {
            descs.push_back(SideShapeFunctionDescriptor(sideIdx, k));
        }
    }
    for (int k = 2; k <= p; k++)
    {
        for (int l = 2; l <= p; l++)
        {
            descs.push_back(InternalShapeFunctionDescriptor(k, l));
        }
    }

    for (const auto& desc : descs)
    {
        m_shapeFunctions[ElementType_Parallelogram].emplace(desc, Polynomial2D(0));
        m_shapeFunctionDerivatives[ElementType_Parallelogram].emplace(desc, std::unordered_map<char, Polynomial2D>{});
    }

    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < descs.size(); i++)
    {
        const auto& desc = descs.at(i);
        Polynomial2D shapeFn = computeQuadShapeFunction(desc);
        Polynomial2D shapeFnDx = diff(shapeFn, 'x');
        Polynomial2D shapeFnDy = diff(shapeFn, 'y');
        m_shapeFunctions[ElementType_Parallelogram].at(desc) = std::move(shapeFn);
        m_shapeFunctionDerivatives[ElementType_Parallelogram].at(desc).emplace('x', std::move(shapeFnDx));
        m_shapeFunctionDerivatives[ElementType_Parallelogram].at(desc).emplace('y', std::move(shapeFnDy));
    }
}

void ShapeFunctionFactory::createTriShapeFunctions(int p)
{
    createLegendrePolynomials(p);
    createShiftedLegendrePolynomials2D(p);
    createRhos(p);

    /* The nodal shape functions are needed for the side and internal shape functions so they must be computed in advance. */
    for (int nodeIdx = 0; nodeIdx < 3; nodeIdx++)
    {
        const NodalShapeFunctionDescriptor desc(nodeIdx);
        Polynomial2D shapeFn = computeTriShapeFunction(desc);
        Polynomial2D shapeFnDx = diff(shapeFn, 'x');
        Polynomial2D shapeFnDy = diff(shapeFn, 'y');
        m_shapeFunctions[ElementType_Triangle].emplace(desc, std::move(shapeFn));
        m_shapeFunctionDerivatives[ElementType_Triangle].emplace(desc, std::unordered_map<char, Polynomial2D>{});
        m_shapeFunctionDerivatives[ElementType_Triangle].at(desc).emplace('x', std::move(shapeFnDx));
        m_shapeFunctionDerivatives[ElementType_Triangle].at(desc).emplace('y', std::move(shapeFnDy));
    }

    std::vector<ShapeFunctionDescriptor> descs;
    for (int sideIdx = 0; sideIdx < 3; sideIdx++)
    {
        for (int k = 2; k <= p; k++)
        {
            descs.push_back(SideShapeFunctionDescriptor(sideIdx, k));
        }
    }
    for (int k = 0; k <= p-2; k++)
    {
        for (int l = 0; l <= p-2; l++)
        {
            descs.push_back(InternalShapeFunctionDescriptor(k, l));
        }
    }

    for (const auto& desc : descs)
    {
        m_shapeFunctions[ElementType_Triangle].emplace(desc, Polynomial2D(0));
        m_shapeFunctionDerivatives[ElementType_Triangle].emplace(desc, std::unordered_map<char, Polynomial2D>{});
    }

    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < descs.size(); i++)
    {
        const auto& desc = descs.at(i);
        Polynomial2D shapeFn = computeTriShapeFunction(desc);
        Polynomial2D shapeFnDx = diff(shapeFn, 'x');
        Polynomial2D shapeFnDy = diff(shapeFn, 'y');
        m_shapeFunctions[ElementType_Triangle].at(desc) = std::move(shapeFn);
        m_shapeFunctionDerivatives[ElementType_Triangle].at(desc).emplace('x', std::move(shapeFnDx));
        m_shapeFunctionDerivatives[ElementType_Triangle].at(desc).emplace('y', std::move(shapeFnDy));
    }
}

void ShapeFunctionFactory::createLegendrePolynomials(int p)
{
    for (int n = 0; n <= p; n++)
    {
        m_legendrePolynomials.emplace(n, computeLegendrePolynomial(n));
    }
}

void ShapeFunctionFactory::createShiftedLegendrePolynomials2D(int p)
{
    for (int n = 0; n <= p - 2; n++)
    {
        m_shiftedLegendrePolynomials2D.emplace(n, std::unordered_map<char, Polynomial2D>{});
    }

    #pragma omp parallel for schedule(static,1)
    for (int n = 0; n <= p - 2; n++)
    {
        m_shiftedLegendrePolynomials2D.at(n).emplace('x', computeShiftedLegendrePolynomial2D(n, 'x'));
        m_shiftedLegendrePolynomials2D.at(n).emplace('y', computeShiftedLegendrePolynomial2D(n, 'y'));
    }
}

void ShapeFunctionFactory::createPhis(int p)
{
    for (int k = 2; k <= p; k++)
    {
        m_phis.emplace(k, Polynomial1D(0));
        m_phis2D.emplace(k, std::unordered_map<char, Polynomial2D>{});
    }

    #pragma omp parallel for schedule(static,1)
    for (int k = 2; k <= p; k++)
    {
        Polynomial1D phi = computePhi(k);
        Polynomial2D phi_x = compose(phi, Polynomial2D("x"));
        Polynomial2D phi_y = compose(phi, Polynomial2D("y"));
        m_phis.at(k) = std::move(phi);
        m_phis2D.at(k).emplace('x', std::move(phi_x));
        m_phis2D.at(k).emplace('y', std::move(phi_y));
    }
}

void ShapeFunctionFactory::createRhos(int p)
{
    for (int k = 2; k <= p; k++)
    {
        m_rhos.emplace(k, Polynomial1D(0));
    }

    #pragma omp parallel for schedule(static,1)
    for (int k = 2; k <= p; k++)
    {
        m_rhos.at(k) = computeRho(k);
    }
}

Polynomial2D ShapeFunctionFactory::computeQuadShapeFunction(const ShapeFunctionDescriptor& d) const
{
    return std::visit([this](const auto& arg) -> Polynomial2D { return this->computeQuadShapeFunction(arg); }, d);
}

Polynomial2D ShapeFunctionFactory::computeQuadShapeFunction(const NodalShapeFunctionDescriptor& d) const
{
    assert(d.nodeIdx < 4);
    if (d.nodeIdx == 0)
    {
        return mpq_class("1/4") * Polynomial2D("1-x") * Polynomial2D("1-y");
    }
    else if (d.nodeIdx == 1)
    {
        return mpq_class("1/4") * Polynomial2D("1+x") * Polynomial2D("1-y");
    }
    else if (d.nodeIdx == 2)
    {
        return mpq_class("1/4") * Polynomial2D("1+x") * Polynomial2D("1+y");
    }
    else
    {
        return mpq_class("1/4") * Polynomial2D("1-x") * Polynomial2D("1+y");
    }
}

Polynomial2D ShapeFunctionFactory::computeQuadShapeFunction(const SideShapeFunctionDescriptor& d) const
{
    assert(d.sideIdx < 4);
    assert(d.k >= 2);
    if (d.sideIdx == 0)
    {
        return (mpq_class("1/2") * Polynomial2D("1-y")) * m_phis2D.at(d.k).at('x');
    }
    else if (d.sideIdx == 1)
    {
        return (mpq_class("1/2") * Polynomial2D("1+x")) * m_phis2D.at(d.k).at('y');
    }
    else if (d.sideIdx == 2)
    {
        return (mpq_class("1/2") * Polynomial2D("1+y")) * compose(m_phis.at(d.k), Polynomial2D("-x"));
    }
    else
    {
        return (mpq_class("1/2") * Polynomial2D("1-x")) * compose(m_phis.at(d.k), Polynomial2D("-y"));
    }
}

Polynomial2D ShapeFunctionFactory::computeQuadShapeFunction(const InternalShapeFunctionDescriptor& d) const
{
    assert(d.k >= 2);
    assert(d.l >= 2);
    return m_phis2D.at(d.k).at('x') * m_phis2D.at(d.l).at('y');
}

Polynomial2D ShapeFunctionFactory::computeTriShapeFunction(const ShapeFunctionDescriptor& d) const
{
    return std::visit([this](const auto& arg) -> Polynomial2D { return this->computeTriShapeFunction(arg); }, d);
}

Polynomial2D ShapeFunctionFactory::computeTriShapeFunction(const NodalShapeFunctionDescriptor& d) const
{
    assert(d.nodeIdx < 3);
    if (d.nodeIdx == 0)
    {
        return Polynomial2D("1-x-y");
    }
    else if (d.nodeIdx == 1)
    {
        return Polynomial2D("x");
    }
    else
    {
        return Polynomial2D("y");
    }
}

Polynomial2D ShapeFunctionFactory::computeTriShapeFunction(const SideShapeFunctionDescriptor& d) const
{
    assert(d.sideIdx < 3);
    assert(d.k >= 2);
    const Polynomial1D& rho = m_rhos.at(d.k);
    const Polynomial2D& nodal1 = m_shapeFunctions[ElementType_Triangle].at(NodalShapeFunctionDescriptor(d.sideIdx));
    const Polynomial2D& nodal2 = m_shapeFunctions[ElementType_Triangle].at(NodalShapeFunctionDescriptor((d.sideIdx + 1) % 3));
    return (nodal1 * nodal2) * compose(rho, nodal2 - nodal1);
}

Polynomial2D ShapeFunctionFactory::computeTriShapeFunction(const InternalShapeFunctionDescriptor& d) const
{
    const Polynomial2D& nodal1 = m_shapeFunctions[ElementType_Triangle].at(NodalShapeFunctionDescriptor(0));
    const Polynomial2D& nodal2 = m_shapeFunctions[ElementType_Triangle].at(NodalShapeFunctionDescriptor(1));
    const Polynomial2D& nodal3 = m_shapeFunctions[ElementType_Triangle].at(NodalShapeFunctionDescriptor(2));
    const Polynomial2D& P_k = m_shiftedLegendrePolynomials2D.at(d.k).at('x');
    const Polynomial2D& P_l = m_shiftedLegendrePolynomials2D.at(d.l).at('y');
    return ((nodal1 * nodal2 * nodal3) * P_k) * P_l;
}

Polynomial1D ShapeFunctionFactory::computeLegendrePolynomial(uint32_t n) const
{
    if (n == 0)
    {
        return Polynomial1D(1);
    }
    else if (n == 1)
    {
        return Polynomial1D("t");
    }
    else
    {
        const mpq_class c1 = mpq_class(2*n - 1) / mpq_class(n);
        const mpq_class c2 = mpq_class(n - 1) / mpq_class(n);
        const Polynomial1D& PmMinus1 = m_legendrePolynomials.at(n-1);
        const Polynomial1D& PmMinus2 = m_legendrePolynomials.at(n-2);
        return (c1 * Polynomial1D("t")) * PmMinus1 - c2 * PmMinus2;
    }
}

Polynomial2D ShapeFunctionFactory::computeShiftedLegendrePolynomial2D(uint32_t n, char variable) const
{
    assert(variable == 'x' || variable == 'y');
    const Polynomial2D g("2" + std::string{variable} + "-1");
    return compose(m_legendrePolynomials.at(n), g);
}

Polynomial1D ShapeFunctionFactory::computePhi(uint32_t k) const
{
    assert(k >= 2);
    const mpq_class C = mpq_class((k-1) * k) / mpq_class(2*k - 1);
    return C * (m_legendrePolynomials.at(k) - m_legendrePolynomials.at(k-2));
}

Polynomial1D ShapeFunctionFactory::computeRho(uint32_t k) const
{
    assert(k >= 2);
    return -4 * diff(m_legendrePolynomials.at(k-1));
}
} // namespace fem
