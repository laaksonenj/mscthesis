#include "fem/basis/TrialFunction.hpp"

#include "fem/basis/BasisFunctionFactory.hpp"
#include "fem/basis/BasisFunctionIndexer.hpp"
#include "fem/basis/ShapeFunctionIndexer.hpp"

namespace fem
{
mpq_class evaluateTrialFunction(const VectorXmpq& coefficients, const FemContext& ctx, const Vector2mpq& x, const ShapeFunctionEvaluator& shapeFunctionEvaluator)
{
    mpq_class res = 0;
    const BasisFunctionIndexer basisFunctionIndexer(ctx);
    const ShapeFunctionIndexer shapeFunctionIndexer(ctx.p, ctx.polynomialSpaceType);
    const Mesh& mesh = *ctx.mesh;
    const Mesh::ElementIndex elementIdx = getIndexOfElementContainingPoint(mesh, x);
    const Element& element = mesh.getElement(elementIdx);
    const ElementType elementType = element.getElementType();
    const AffineMap Finv = element.getReferenceElementMap().inverse();
    for (int shapeFnIdx = 0; shapeFnIdx < basisFunctionIndexer.getNumOfShapeFunctions(elementIdx); shapeFnIdx++)
    {
        const uint32_t basisFnIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFnIdx);
        const auto desc = shapeFunctionIndexer.getShapeFunctionDescriptor(elementType, shapeFnIdx);
        mpq_class value = shapeFunctionEvaluator.evaluate(elementType, desc, Finv(x));
        if (const auto* descp = std::get_if<SideShapeFunctionDescriptor>(&desc))
        {
            const auto adjacentElementIdx = mesh.getIndexOfAdjacentElement(elementIdx, descp->sideIdx);
            if (adjacentElementIdx.has_value() && elementIdx < adjacentElementIdx.value() && descp->k % 2 != 0)
            {
                value *= -1;
            }
        }
        res += coefficients(basisFnIdx) * value;
    }
    return res;
}

mpq_class integrateTrialFunction(const VectorXmpq& coefficients, const FemContext& ctx, const ShapeFunctionFactory& shapeFunctionFactory)
{
    mpq_class res = 0;
    const BasisFunctionIndexer basisFunctionIndexer(ctx);
    const ShapeFunctionIndexer shapeFunctionIndexer(ctx.p, ctx.polynomialSpaceType);
    const Mesh& mesh = *ctx.mesh;
    for (int elementIdx = 0; elementIdx < mesh.getNumOfElements(); elementIdx++)
    {
        const Element& element = mesh.getElement(elementIdx);
        const ElementType elementType = element.getElementType();
        const mpq_class detA = element.getReferenceElementMap().A.determinant();
        for (int shapeFnIdx = 0; shapeFnIdx < shapeFunctionIndexer.getNumOfShapeFunctions(elementType); shapeFnIdx++)
        {
            const uint32_t basisFnIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFnIdx);
            const auto desc = shapeFunctionIndexer.getShapeFunctionDescriptor(elementType, shapeFnIdx);
            const Polynomial2D& shapeFn = shapeFunctionFactory.getShapeFunction(elementType, desc);
            mpq_class integral = integrateOverReferenceElement(shapeFn, elementType);
            if (const auto* descp = std::get_if<SideShapeFunctionDescriptor>(&desc))
            {
                const auto adjacentElementIdx = mesh.getIndexOfAdjacentElement(elementIdx, descp->sideIdx);
                if (adjacentElementIdx.has_value() && elementIdx < adjacentElementIdx.value() && descp->k % 2 != 0)
                {
                    integral *= -1; // is this even needed? integral may always be zero
                }
            }
            res += coefficients(basisFnIdx) * detA * integral;
        }
    }
    return res;
}

void normalizeTrialFunction(VectorXmpq& coefficients, const FemContext& ctx, const ShapeFunctionFactory& shapeFunctionFactory)
{
    const Mesh& mesh = *ctx.mesh;
    const mpq_class areaOfMesh = calculateMeshArea(mesh);
    const mpq_class normalizationConst = integrateTrialFunction(coefficients, ctx, shapeFunctionFactory) / areaOfMesh;
    for (int nodeIdx = 0; nodeIdx < mesh.getNumOfNodes(); nodeIdx++)
    {
        coefficients(nodeIdx) -= normalizationConst;
    }
}

mpq_class evaluateTrialFunction(const VectorXmpq& coefficients, const FemContext& ctx, const Vector2mpq& x)
{
    mpq_class res = 0;
    const BasisFunctionFactory basisFunctionFactory(ctx);
    const BasisFunctionIndexer basisFunctionIndexer(ctx);
    const Mesh& mesh = *ctx.mesh;
    const Mesh::ElementIndex elementIdx = getIndexOfElementContainingPoint(mesh, x);
    const Element& element = mesh.getElement(elementIdx);
    const AffineMap Finv = element.getReferenceElementMap().inverse();
    for (int shapeFunctionIdx = 0; shapeFunctionIdx < basisFunctionIndexer.getNumOfShapeFunctions(elementIdx); shapeFunctionIdx++)
    {
        const uint32_t basisFunctionIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
        const Polynomial2D shapeFn = basisFunctionFactory.getShapeFunction(elementIdx, shapeFunctionIdx);
        res += coefficients(basisFunctionIdx) * shapeFn(Finv(x));
    }
    return res;
}

mpq_class integrateTrialFunction(const VectorXmpq& coefficients, const FemContext& ctx)
{
    mpq_class res = 0;
    const BasisFunctionFactory basisFunctionFactory(ctx);
    const BasisFunctionIndexer basisFunctionIndexer(ctx);
    const Mesh& mesh = *ctx.mesh;
    for (int elementIdx = 0; elementIdx < mesh.getNumOfElements(); elementIdx++)
    {
        const Element& element = mesh.getElement(elementIdx);
        const AffineMap F = element.getReferenceElementMap();
        const mpq_class detA = F.A.determinant();
        for (int shapeFunctionIdx = 0; shapeFunctionIdx < basisFunctionIndexer.getNumOfShapeFunctions(elementIdx); shapeFunctionIdx++)
        {
            const uint32_t basisFunctionIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
            const Polynomial2D shapeFn = basisFunctionFactory.getShapeFunction(elementIdx, shapeFunctionIdx);
            res += coefficients(basisFunctionIdx) * detA * integrateOverReferenceElement(shapeFn, element.getElementType());
        }
    }
    return res;
}
} // namespace fem
