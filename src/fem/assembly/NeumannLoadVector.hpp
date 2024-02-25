#pragma once

#include "fem/basis/FemContext.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/math/Function.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
VectorXmpq assembleNeumannLoadVector(const FemContext& ctx,
                                     const BivariateFunction& g,
                                     Mesh::ElementIndex elementIdx,
                                     Mesh::SideIndex localSideIdx,
                                     const ShapeFunctionFactory& shapeFunctionFactory);
                                     
VectorXmpq assembleNeumannLoadVector(const FemContext& ctx,
                                     const GradientFunction& grad_u,
                                     const ShapeFunctionFactory& shapeFunctionFactory);

/* For testing purposes */
VectorXmpq assembleNeumannLoadVector(const FemContext& ctx,
                                     const BivariateFunction& g,
                                     Mesh::ElementIndex elementIdx,
                                     Mesh::SideIndex localSideIdx);
} // namespace fem
