#include "fem/basis/ShapeFunctionFactory.hpp"

#include <cassert>
#include <vector>

namespace fem
{
ShapeFunctionFactory::ShapeFunctionFactory()
    : m_numOfPrecomputedQuad(0)
    , m_numOfPrecomputedTri(0)
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

void ShapeFunctionFactory::precomputeShapeFunctions(ElementType elementType, uint32_t p)
{
    assert(p >= 1);
    if (elementType == ElementType_Parallelogram)
    {
        if (p > m_numOfPrecomputedQuad)
        {
            precomputeQuadShapeFunctions(p);
            m_numOfPrecomputedQuad = p;
        }
    }
    else
    {
        if (p > m_numOfPrecomputedTri)
        {
            precomputeTriShapeFunctions(p);
            m_numOfPrecomputedTri = p;
        }
    }
}

const Polynomial2D& ShapeFunctionFactory::getQuadShapeFunction(const ShapeFunctionDescriptor& d) const
{
    if (!m_quadShapeFunctions.contains(d))
    {
        m_quadShapeFunctions.emplace(d, computeQuadShapeFunction(d));
    }
    return m_quadShapeFunctions.at(d);
}

const Polynomial2D& ShapeFunctionFactory::getTriShapeFunction(const ShapeFunctionDescriptor& d) const
{
    if (!m_triShapeFunctions.contains(d))
    {
        m_triShapeFunctions.emplace(d, computeTriShapeFunction(d));
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
        m_quadShapeFunctionDerivatives.at(d).emplace(variable, computeQuadShapeFunctionDerivative(d, variable));
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
        m_triShapeFunctionDerivatives.at(d).emplace(variable, computeTriShapeFunctionDerivative(d, variable));
    }
    return m_triShapeFunctionDerivatives.at(d).at(variable);
}

const Polynomial1D& ShapeFunctionFactory::getLegendrePolynomial(uint32_t n) const
{
    if (!m_legendrePolynomials.contains(n))
    {
        m_legendrePolynomials.emplace(n, computeLegendrePolynomial(n));
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
        m_shiftedLegendrePolynomials2D.at(n).emplace(variable, computeShiftedLegendrePolynomial2D(n, variable));
    }
    return m_shiftedLegendrePolynomials2D.at(n).at(variable);
}
    
const Polynomial1D& ShapeFunctionFactory::getPhi(uint32_t k) const
{
    assert(k >= 2);
    if (!m_phis.contains(k))
    {
        m_phis.emplace(k, computePhi(k));
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
        m_phis2D.at(k).emplace(variable, computePhi2D(k, variable));
    }
    return m_phis2D.at(k).at(variable);
}

const Polynomial1D& ShapeFunctionFactory::getRho(uint32_t k) const
{
    assert(k >= 2);
    if (!m_rhos.contains(k))
    {
        m_rhos.emplace(k, computeRho(k));
    }
    return m_rhos.at(k);
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
        return (mpq_class("1/2") * Polynomial2D("1-y")) * getPhi2D(d.k, 'x');
    }
    else if (d.sideIdx == 1)
    {
        return (mpq_class("1/2") * Polynomial2D("1+x")) * getPhi2D(d.k, 'y');
    }
    else if (d.sideIdx == 2)
    {
        return (mpq_class("1/2") * Polynomial2D("1+y")) * compose(getPhi(d.k), Polynomial2D("-x"));
    }
    else
    {
        return (mpq_class("1/2") * Polynomial2D("1-x")) * compose(getPhi(d.k), Polynomial2D("-y"));
    }
}

Polynomial2D ShapeFunctionFactory::computeQuadShapeFunction(const InternalShapeFunctionDescriptor& d) const
{
    assert(d.k >= 2);
    assert(d.l >= 2);
    return getPhi2D(d.k, 'x') * getPhi2D(d.l, 'y');
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
    const Polynomial1D& rho = getRho(d.k);
    const Polynomial2D& nodal1 = getTriShapeFunction(NodalShapeFunctionDescriptor(d.sideIdx));
    const Polynomial2D& nodal2 = getTriShapeFunction(NodalShapeFunctionDescriptor((d.sideIdx + 1) % 3));
    return (nodal1 * nodal2) * compose(rho, nodal2 - nodal1);
}

Polynomial2D ShapeFunctionFactory::computeTriShapeFunction(const InternalShapeFunctionDescriptor& d) const
{
    const Polynomial2D& nodal1 = getTriShapeFunction(NodalShapeFunctionDescriptor(0));
    const Polynomial2D& nodal2 = getTriShapeFunction(NodalShapeFunctionDescriptor(1));
    const Polynomial2D& nodal3 = getTriShapeFunction(NodalShapeFunctionDescriptor(2));
    const Polynomial2D& P_k = getShiftedLegendrePolynomial2D(d.k, 'x');
    const Polynomial2D& P_l = getShiftedLegendrePolynomial2D(d.l, 'y');
    return ((nodal1 * nodal2 * nodal3) * P_k) * P_l;
}

Polynomial2D ShapeFunctionFactory::computeQuadShapeFunctionDerivative(const ShapeFunctionDescriptor& d, char variable) const
{
    assert(variable == 'x' || variable == 'y');
    return diff(getQuadShapeFunction(d), variable);
}

Polynomial2D ShapeFunctionFactory::computeTriShapeFunctionDerivative(const ShapeFunctionDescriptor& d, char variable) const
{
    assert(variable == 'x' || variable == 'y');
    return diff(getTriShapeFunction(d), variable);
}

Polynomial1D ShapeFunctionFactory::computeLegendrePolynomial(uint32_t n) const
{
    const mpq_class c1 = mpq_class(2*n - 1) / mpq_class(n);
    const mpq_class c2 = mpq_class(n - 1) / mpq_class(n);
    const Polynomial1D& PmMinus1 = getLegendrePolynomial(n-1);
    const Polynomial1D& PmMinus2 = getLegendrePolynomial(n-2);
    return (c1 * Polynomial1D("t")) * PmMinus1 - c2 * PmMinus2;
}

Polynomial2D ShapeFunctionFactory::computeShiftedLegendrePolynomial2D(uint32_t n, char variable) const
{
    assert(variable == 'x' || variable == 'y');
    const Polynomial2D g("2" + std::string{variable} + "-1");
    return compose(getLegendrePolynomial(n), g);
}

Polynomial1D ShapeFunctionFactory::computePhi(uint32_t k) const
{
    assert(k >= 2);
    const mpq_class C = mpq_class((k-1) * k) / mpq_class(2*k - 1);
    return C * (getLegendrePolynomial(k) - getLegendrePolynomial(k-2));
}

Polynomial2D ShapeFunctionFactory::computePhi2D(uint32_t k, char variable) const
{
    assert(variable == 'x' || variable == 'y');
    assert(k >= 2);
    return compose(getPhi(k), Polynomial2D(std::string{variable}));
}

Polynomial1D ShapeFunctionFactory::computeRho(uint32_t k) const
{
    assert(k >= 2);
    return -4 * diff(getLegendrePolynomial(k-1));
}

void ShapeFunctionFactory::precomputeQuadShapeFunctions(int p)
{
    precomputeLegendrePolynomials(p);
    precomputePhis(p);

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
        m_quadShapeFunctions.emplace(desc, Polynomial2D(0));
        m_quadShapeFunctionDerivatives.emplace(desc, std::unordered_map<char, Polynomial2D>{});
    }

    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < descs.size(); i++)
    {
        const auto& desc = descs.at(i);
        Polynomial2D shapeFn = computeQuadShapeFunction(desc);
        Polynomial2D shapeFnDx = diff(shapeFn, 'x');
        Polynomial2D shapeFnDy = diff(shapeFn, 'y');
        m_quadShapeFunctions.at(desc) = std::move(shapeFn);
        m_quadShapeFunctionDerivatives.at(desc).emplace('x', std::move(shapeFnDx));
        m_quadShapeFunctionDerivatives.at(desc).emplace('y', std::move(shapeFnDy));
    }
}

void ShapeFunctionFactory::precomputeTriShapeFunctions(int p)
{
    precomputeLegendrePolynomials(p);
    precomputeShiftedLegendrePolynomials2D(p);
    precomputeRhos(p);

    /* The nodal shape functions are needed for the side and internal shape functions so they must be computed in advance. */
    for (int nodeIdx = 0; nodeIdx < 3; nodeIdx++)
    {
        const NodalShapeFunctionDescriptor desc(nodeIdx);
        Polynomial2D shapeFn = computeTriShapeFunction(desc);
        Polynomial2D shapeFnDx = diff(shapeFn, 'x');
        Polynomial2D shapeFnDy = diff(shapeFn, 'y');
        m_triShapeFunctions.emplace(desc, shapeFn);
        m_triShapeFunctionDerivatives.emplace(desc, std::unordered_map<char, Polynomial2D>{});
        m_triShapeFunctionDerivatives.at(desc).emplace('x', shapeFnDx);
        m_triShapeFunctionDerivatives.at(desc).emplace('y', shapeFnDy);
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
        m_triShapeFunctions.emplace(desc, Polynomial2D(0));
        m_triShapeFunctionDerivatives.emplace(desc, std::unordered_map<char, Polynomial2D>{});
    }

    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < descs.size(); i++)
    {
        const auto& desc = descs.at(i);
        Polynomial2D shapeFn = computeTriShapeFunction(desc);
        Polynomial2D shapeFnDx = diff(shapeFn, 'x');
        Polynomial2D shapeFnDy = diff(shapeFn, 'y');
        m_triShapeFunctions.at(desc) = std::move(shapeFn);
        m_triShapeFunctionDerivatives.at(desc).emplace('x', std::move(shapeFnDx));
        m_triShapeFunctionDerivatives.at(desc).emplace('y', std::move(shapeFnDy));
    }
}

void ShapeFunctionFactory::precomputeLegendrePolynomials(int p)
{
    const uint32_t currentNum = m_legendrePolynomials.size();
    for (int n = currentNum; n <= p; n++)
    {
        m_legendrePolynomials.emplace(n, computeLegendrePolynomial(n));
    }
}

void ShapeFunctionFactory::precomputeShiftedLegendrePolynomials2D(int p)
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

void ShapeFunctionFactory::precomputePhis(int p)
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

void ShapeFunctionFactory::precomputeRhos(int p)
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
} // namespace fem
