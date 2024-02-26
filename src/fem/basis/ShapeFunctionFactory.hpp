#pragma once

#include <cstdint>
#include <unordered_map>

#include "fem/basis/ShapeFunctionDescriptor.hpp"
#include "fem/domain/Element.hpp"
#include "fem/math/Polynomial.hpp"

namespace fem
{
class ShapeFunctionFactory
{
public:
    ShapeFunctionFactory() = default;

    void createShapeFunctions(ElementType elementType, uint32_t p);

    const Polynomial2D& getShapeFunction(ElementType elementType, const ShapeFunctionDescriptor& descriptor) const;
    const Polynomial2D& getShapeFunctionDerivative(ElementType elementType, const ShapeFunctionDescriptor& descriptor, char variable) const;

private:

    void createQuadShapeFunctions(int p);
    void createTriShapeFunctions(int p);
    void createLegendrePolynomials(int p);
    void createShiftedLegendrePolynomials2D(int p);
    void createPhis(int p);
    void createRhos(int p);

    Polynomial2D computeQuadShapeFunction(const ShapeFunctionDescriptor& d) const;
    Polynomial2D computeQuadShapeFunction(const NodalShapeFunctionDescriptor& d) const;
    Polynomial2D computeQuadShapeFunction(const SideShapeFunctionDescriptor& d) const;
    Polynomial2D computeQuadShapeFunction(const InternalShapeFunctionDescriptor& d) const;

    Polynomial2D computeTriShapeFunction(const ShapeFunctionDescriptor& d) const;
    Polynomial2D computeTriShapeFunction(const NodalShapeFunctionDescriptor& d) const;
    Polynomial2D computeTriShapeFunction(const SideShapeFunctionDescriptor& d) const;
    Polynomial2D computeTriShapeFunction(const InternalShapeFunctionDescriptor& d) const;

    Polynomial1D computeLegendrePolynomial(uint32_t n) const;
    Polynomial2D computeShiftedLegendrePolynomial2D(uint32_t n, char variable) const;
    Polynomial1D computePhi(uint32_t k) const;
    Polynomial1D computeRho(uint32_t k) const;

private:
    std::unordered_map<ShapeFunctionDescriptor, Polynomial2D> m_shapeFunctions[2];
    std::unordered_map<ShapeFunctionDescriptor, std::unordered_map<char, Polynomial2D>> m_shapeFunctionDerivatives[2];

    std::unordered_map<uint32_t, Polynomial1D> m_legendrePolynomials;
    std::unordered_map<uint32_t, std::unordered_map<char, Polynomial2D>> m_shiftedLegendrePolynomials2D;
    std::unordered_map<uint32_t, Polynomial1D> m_phis;
    std::unordered_map<uint32_t, std::unordered_map<char, Polynomial2D>> m_phis2D;
    std::unordered_map<uint32_t, Polynomial1D> m_rhos;
};
} // namespace fem
