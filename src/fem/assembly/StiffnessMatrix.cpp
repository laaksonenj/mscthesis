#include "fem/assembly/StiffnessMatrix.hpp"

#include <algorithm>
#include <cassert>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <boost/functional/hash.hpp>

#include "fem/basis/BasisFunctionFactory.hpp"
#include "fem/basis/BasisFunctionIndexer.hpp"
#include "fem/basis/ShapeFunctionIndexer.hpp"

namespace fem
{
namespace
{
struct StiffnessMatrixAssembler
{
    using ShapeFunctionDerivativePair = std::tuple<uint32_t, char, uint32_t, char>;
    using ShapeFunctionIntegralCache = std::unordered_map<ShapeFunctionDerivativePair, mpq_class, boost::hash<ShapeFunctionDerivativePair>>;

    const FemContext& ctx;
    const ShapeFunctionFactory& shapeFunctionFactory;
    BasisFunctionIndexer basisFunctionIndexer;
    ShapeFunctionIndexer shapeFunctionIndexer;
    MatrixXmpq stiffnessMatrix;
    ShapeFunctionIntegralCache integralCache[2]; // one for triangles and one for quads

    StiffnessMatrixAssembler(const FemContext& ctx, const ShapeFunctionFactory& shapeFunctionFactory);
    void assemble();
    void precomputeIntegrals(ElementType elementType);
    void assembleElement(Mesh::ElementIndex elementIdx);
    void symmetrize();
};

StiffnessMatrixAssembler::StiffnessMatrixAssembler(const FemContext& ctx, const ShapeFunctionFactory& shapeFunctionFactory)
    : ctx(ctx)
    , shapeFunctionFactory(shapeFunctionFactory)
    , basisFunctionIndexer(BasisFunctionIndexer(ctx))
    , shapeFunctionIndexer(ShapeFunctionIndexer(ctx.p, ctx.polynomialSpaceType))
{
    const uint32_t numOfBasisFunctions = basisFunctionIndexer.getNumOfBasisFunctions();
    stiffnessMatrix = MatrixXmpq(numOfBasisFunctions, numOfBasisFunctions);
}

void StiffnessMatrixAssembler::assemble()
{
    const Mesh& mesh = *ctx.mesh;
    if (mesh.containsQuadrilateral())
    {
        precomputeIntegrals(ElementType_Parallelogram);
    }
    if (mesh.containsTriangle())
    {
        precomputeIntegrals(ElementType_Triangle);
    }
    for (int elementIdx = 0; elementIdx < ctx.mesh->getNumOfElements(); elementIdx++)
    {
        assembleElement(elementIdx);
    }
    symmetrize();
}

void StiffnessMatrixAssembler::assembleElement(Mesh::ElementIndex elementIdx)
{
    const Mesh& mesh = *ctx.mesh;
    const Element& element = mesh.getElement(elementIdx);
    const ElementType elementType = element.getElementType();
    const AffineMap F = element.getReferenceElementMap();
    const Matrix2mpq& A = F.A;
    const Matrix2mpq Ainv = A.inverse();
    const Matrix2mpq AinvT = Ainv.transpose();
    const Matrix2mpq M = Ainv * AinvT;
    const mpq_class detA = A.determinant();
    const uint32_t numOfShapeFunctions = shapeFunctionIndexer.getNumOfShapeFunctions(elementType);
    for (int shapeFnIdx1 = 0; shapeFnIdx1 < numOfShapeFunctions; shapeFnIdx1++)
    {
        const auto desc1 = shapeFunctionIndexer.getShapeFunctionDescriptor(elementType, shapeFnIdx1);
        const uint32_t basisFnIdx1 = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFnIdx1);
        for (int shapeFnIdx2 = shapeFnIdx1; shapeFnIdx2 < numOfShapeFunctions; shapeFnIdx2++)
        {
            const auto desc2 = shapeFunctionIndexer.getShapeFunctionDescriptor(elementType, shapeFnIdx2);
            const uint32_t basisFnIdx2 = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFnIdx2);

            const auto key_xx = std::make_tuple(shapeFnIdx1, 'x', shapeFnIdx2, 'x');
            const auto key_xy = std::make_tuple(shapeFnIdx1, 'x', shapeFnIdx2, 'y');
            const auto key_yx = std::make_tuple(shapeFnIdx1, 'y', shapeFnIdx2, 'x');
            const auto key_yy = std::make_tuple(shapeFnIdx1, 'y', shapeFnIdx2, 'y');

            mpq_class integral = 0;
            integral += M(0,0) * integralCache[elementType].at(key_xx);
            integral += M(0,1) * integralCache[elementType].at(key_xy);
            integral += M(1,0) * integralCache[elementType].at(key_yx);
            integral += M(1,1) * integralCache[elementType].at(key_yy);
            integral *= detA;

            if (const auto* descp1 = std::get_if<SideShapeFunctionDescriptor>(&desc1))
            {
                const auto adjacentElementIdx = mesh.getIndexOfAdjacentElement(elementIdx, descp1->sideIdx);
                if (adjacentElementIdx.has_value() && elementIdx < adjacentElementIdx.value() && descp1->k % 2 != 0)
                {
                    integral *= -1;
                }
            }
            if (const auto* descp2 = std::get_if<SideShapeFunctionDescriptor>(&desc2))
            {
                const auto adjacentElementIdx = mesh.getIndexOfAdjacentElement(elementIdx, descp2->sideIdx);
                if (adjacentElementIdx.has_value() && elementIdx < adjacentElementIdx.value() && descp2->k % 2 != 0)
                {
                    integral *= -1;
                }
            }

            const uint32_t i = std::min(basisFnIdx1, basisFnIdx2);
            const uint32_t j = std::max(basisFnIdx1, basisFnIdx2);
            stiffnessMatrix(i,j) += std::move(integral);
        }
    }
}

void StiffnessMatrixAssembler::symmetrize()
{
    const uint32_t dim = stiffnessMatrix.rows();
    for (int j = 0; j < dim; j++)
    {
        for (int i = j + 1; i < dim; i++)
        {
            stiffnessMatrix(i,j) = stiffnessMatrix(j,i);
        }
    }
}

void StiffnessMatrixAssembler::precomputeIntegrals(ElementType elementType)
{
    auto& cache = integralCache[elementType];
    const uint32_t numOfShapeFunctions = shapeFunctionIndexer.getNumOfShapeFunctions(elementType);
    std::vector<ShapeFunctionDerivativePair> derivativePairs;
    for (int shapeFnIdx1 = 0; shapeFnIdx1 < numOfShapeFunctions; shapeFnIdx1++)
    {
        for (int shapeFnIdx2 = shapeFnIdx1; shapeFnIdx2 < numOfShapeFunctions; shapeFnIdx2++)
        {
            const auto key_xx = std::make_tuple(shapeFnIdx1, 'x', shapeFnIdx2, 'x');
            const auto key_xy = std::make_tuple(shapeFnIdx1, 'x', shapeFnIdx2, 'y');
            const auto key_yx = std::make_tuple(shapeFnIdx1, 'y', shapeFnIdx2, 'x');
            const auto key_yy = std::make_tuple(shapeFnIdx1, 'y', shapeFnIdx2, 'y');
            cache.emplace(key_xx, 0);
            cache.emplace(key_xy, 0);
            cache.emplace(key_yx, 0);
            cache.emplace(key_yy, 0);
            derivativePairs.push_back(key_xx);
            derivativePairs.push_back(key_xy);
            derivativePairs.push_back(key_yx);
            derivativePairs.push_back(key_yy);
        }
    }

    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < derivativePairs.size(); i++)
    {
        const auto& derivativePair = derivativePairs[i];
        const uint32_t shapeFnIdx1 = std::get<0>(derivativePair);
        const char var1 = std::get<1>(derivativePair);
        const uint32_t shapeFnIdx2 = std::get<2>(derivativePair);
        const char var2 = std::get<3>(derivativePair);
        const auto desc1 = shapeFunctionIndexer.getShapeFunctionDescriptor(elementType, shapeFnIdx1);
        const auto desc2 = shapeFunctionIndexer.getShapeFunctionDescriptor(elementType, shapeFnIdx2);
        const Polynomial2D& shapeFn1D = shapeFunctionFactory.getShapeFunctionDerivative(elementType, desc1, var1);
        const Polynomial2D& shapeFn2D = shapeFunctionFactory.getShapeFunctionDerivative(elementType, desc2, var2);
        cache.at(derivativePair) = integrateOverReferenceElement(shapeFn1D * shapeFn2D, elementType);
    }
}
} // namespace

MatrixXmpq assembleStiffnessMatrix(const FemContext& ctx, const ShapeFunctionFactory& shapeFunctionFactory)
{
    StiffnessMatrixAssembler assembler(ctx, shapeFunctionFactory);
    assembler.assemble();
    return std::move(assembler.stiffnessMatrix);
}

MatrixXmpq assembleStiffnessMatrix(const FemContext& ctx)
{
    const BasisFunctionFactory basisFunctionFactory(ctx);
    const BasisFunctionIndexer basisFunctionIndexer(ctx);
    const uint32_t numOfBasisFunctions = basisFunctionIndexer.getNumOfBasisFunctions();
    MatrixXmpq res(numOfBasisFunctions, numOfBasisFunctions);
    const Mesh& mesh = *ctx.mesh;
    for (int elementIdx = 0; elementIdx < mesh.getNumOfElements(); elementIdx++)
    {
        const Element& element = mesh.getElement(elementIdx);
        const AffineMap F = element.getReferenceElementMap();
        const Matrix2mpq& A = F.A;
        const Matrix2mpq Ainv = A.inverse();
        const Matrix2mpq AinvT = Ainv.transpose();
        const Matrix2mpq M = Ainv * AinvT;
        const mpq_class detA = abs(A.determinant());
        const uint32_t numOfShapeFunctions = basisFunctionIndexer.getNumOfShapeFunctions(elementIdx);
        for (int shapeFunctionIdx1 = 0; shapeFunctionIdx1 < numOfShapeFunctions; shapeFunctionIdx1++)
        {
            for (int shapeFunctionIdx2 = 0; shapeFunctionIdx2 < numOfShapeFunctions; shapeFunctionIdx2++)
            {
                const Polynomial2D shapeFunction1Dx = basisFunctionFactory.getShapeFunctionDerivative(elementIdx, shapeFunctionIdx1, 'x');
                const Polynomial2D shapeFunction1Dy = basisFunctionFactory.getShapeFunctionDerivative(elementIdx, shapeFunctionIdx1, 'y');
                const Polynomial2D shapeFunction2Dx = basisFunctionFactory.getShapeFunctionDerivative(elementIdx, shapeFunctionIdx2, 'x');
                const Polynomial2D shapeFunction2Dy = basisFunctionFactory.getShapeFunctionDerivative(elementIdx, shapeFunctionIdx2, 'y');
                mpq_class integral = 0;
                integral += M(0,0) * integrateOverReferenceElement(shapeFunction1Dx * shapeFunction2Dx, element.getElementType());
                integral += M(0,1) * integrateOverReferenceElement(shapeFunction1Dx * shapeFunction2Dy, element.getElementType());
                integral += M(1,0) * integrateOverReferenceElement(shapeFunction1Dy * shapeFunction2Dx, element.getElementType());
                integral += M(1,1) * integrateOverReferenceElement(shapeFunction1Dy * shapeFunction2Dy, element.getElementType());
                integral *= detA;
                const uint32_t basisFunctionIdx1 = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx1);
                const uint32_t basisFunctionIdx2 = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx2);
                res(basisFunctionIdx1, basisFunctionIdx2) += integral;
            }
        }
    }
    return res;
}

MatrixXmpq extractSubStiffnessMatrix(const MatrixXmpq& stiffnessMatrix, const FemContext& ctx, uint32_t p)
{
    assert(p >= 1 && p <= ctx.p);
    const BasisFunctionIndexer superBasisFunctionIndexer(ctx);
    const BasisFunctionIndexer subBasisFunctionIndexer(FemContext(ctx.mesh, p, ctx.polynomialSpaceType));
    const uint32_t numOfBasisFunctions = subBasisFunctionIndexer.getNumOfBasisFunctions();
    MatrixXmpq res(numOfBasisFunctions, numOfBasisFunctions);
    for (int i = 0; i < numOfBasisFunctions; i++)
    {
        const BasisFunctionDescriptor desc_i = subBasisFunctionIndexer.getBasisFunctionDescriptor(i);
        const int ii = superBasisFunctionIndexer.getBasisFunctionIndex(desc_i);
        for (int j = 0; j < numOfBasisFunctions; j++)
        {
            const BasisFunctionDescriptor desc_j = subBasisFunctionIndexer.getBasisFunctionDescriptor(j);
            const int jj = superBasisFunctionIndexer.getBasisFunctionIndex(desc_j);
            res(i,j) = stiffnessMatrix(ii,jj);
        }
    }
    return res;
}
} // namespace fem
