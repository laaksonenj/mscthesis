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
    ShapeFunctionFactory();

    const Polynomial2D& getShapeFunction(ElementType elementType, const ShapeFunctionDescriptor& descriptor) const;
    const Polynomial2D& getShapeFunctionDerivative(ElementType elementType, const ShapeFunctionDescriptor& descriptor, char variable) const;

private:
    const Polynomial2D& getQuadShapeFunction(const ShapeFunctionDescriptor& d) const;
    const Polynomial2D& getQuadShapeFunction(const NodalShapeFunctionDescriptor& d) const;
    const Polynomial2D& getQuadShapeFunction(const SideShapeFunctionDescriptor& d) const;
    const Polynomial2D& getQuadShapeFunction(const InternalShapeFunctionDescriptor& d) const;

    const Polynomial2D& getTriShapeFunction(const ShapeFunctionDescriptor& d) const;
    const Polynomial2D& getTriShapeFunction(const NodalShapeFunctionDescriptor& d) const;
    const Polynomial2D& getTriShapeFunction(const SideShapeFunctionDescriptor& d) const;
    const Polynomial2D& getTriShapeFunction(const InternalShapeFunctionDescriptor& d) const;

    const Polynomial2D& getQuadShapeFunctionDerivative(const ShapeFunctionDescriptor& d, char variable) const;
    const Polynomial2D& getTriShapeFunctionDerivative(const ShapeFunctionDescriptor& d, char variable) const;

    const Polynomial1D& getLegendrePolynomial(uint32_t n) const;
    const Polynomial2D& getShiftedLegendrePolynomial2D(uint32_t n, char variable) const;
    const Polynomial1D& getPhi(uint32_t k) const;
    const Polynomial2D& getPhi2D(uint32_t k, char variable) const;
    const Polynomial1D& getRho(uint32_t k) const;

private:
    mutable std::unordered_map<ShapeFunctionDescriptor, Polynomial2D> m_quadShapeFunctions;
    mutable std::unordered_map<ShapeFunctionDescriptor, Polynomial2D> m_triShapeFunctions;
    mutable std::unordered_map<ShapeFunctionDescriptor, std::unordered_map<char, Polynomial2D>> m_quadShapeFunctionDerivatives;
    mutable std::unordered_map<ShapeFunctionDescriptor, std::unordered_map<char, Polynomial2D>> m_triShapeFunctionDerivatives;

    mutable std::unordered_map<uint32_t, Polynomial1D> m_legendrePolynomials;
    mutable std::unordered_map<uint32_t, std::unordered_map<char, Polynomial2D>> m_shiftedLegendrePolynomials2D;
    mutable std::unordered_map<uint32_t, Polynomial1D> m_phis;
    mutable std::unordered_map<uint32_t, std::unordered_map<char, Polynomial2D>> m_phis2D;
    mutable std::unordered_map<uint32_t, Polynomial1D> m_rhos;
};
} // namespace fem
