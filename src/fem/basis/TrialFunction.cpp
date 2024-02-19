#include "fem/basis/TrialFunction.hpp"

namespace fem
{
mpq_class evaluateTrialFunction(const VectorXmpq& coefficients, const BasisFunctionIndexer& basisFunctionIndexer, const Vector2mpq& p, const ShapeFunctionFactory& shapeFunctionFactory)
{
    return evaluateTrialFunction(coefficients, basisFunctionIndexer, p);
}

mpq_class integrateTrialFunction(const VectorXmpq& coefficients, const BasisFunctionIndexer& basisFunctionIndexer, const ShapeFunctionFactory& shapeFunctionFactory)
{
    return integrateTrialFunction(coefficients, basisFunctionIndexer);
}

void normalizeTrialFunction(VectorXmpq& coefficients, const BasisFunctionIndexer& basisFunctionIndexer, const ShapeFunctionFactory& shapeFunctionFactory)
{
    const Mesh& mesh = basisFunctionIndexer.getMesh();
    const mpq_class areaOfMesh = calculateMeshArea(mesh);
    const mpq_class normalizationConst = integrateTrialFunction(coefficients, basisFunctionIndexer, shapeFunctionFactory) / areaOfMesh;
    for (int nodeIdx = 0; nodeIdx < mesh.getNumOfNodes(); nodeIdx++)
    {
        coefficients(nodeIdx) -= normalizationConst;
    }
}

mpq_class evaluateTrialFunction(const VectorXmpq& coefficients, const BasisFunctionIndexer& basisFunctionIndexer, const Vector2mpq& p)
{
    mpq_class res = 0;
    const Mesh& mesh = basisFunctionIndexer.getMesh();
    const Mesh::ElementIndex elementIdx = getIndexOfElementContainingPoint(mesh, p);
    const Element& element = mesh.getElement(elementIdx);
    const AffineMap Finv = element.getReferenceElementMap().inverse();
    for (int shapeFunctionIdx = 0; shapeFunctionIdx < basisFunctionIndexer.getNumOfShapeFunctions(elementIdx); shapeFunctionIdx++)
    {
        const uint32_t basisFunctionIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
        const Polynomial2D shapeFn = basisFunctionIndexer.getShapeFunction(elementIdx, shapeFunctionIdx);
        res += coefficients(basisFunctionIdx) * shapeFn(Finv(p));
    }
    return res;
}

mpq_class integrateTrialFunction(const VectorXmpq& coefficients, const BasisFunctionIndexer& basisFunctionIndexer)
{
    mpq_class res = 0;
    const Mesh& mesh = basisFunctionIndexer.getMesh();
    for (int elementIdx = 0; elementIdx < mesh.getNumOfElements(); elementIdx++)
    {
        const Element& element = mesh.getElement(elementIdx);
        const AffineMap F = element.getReferenceElementMap();
        const mpq_class detA = F.A.determinant();
        for (int shapeFunctionIdx = 0; shapeFunctionIdx < basisFunctionIndexer.getNumOfShapeFunctions(elementIdx); shapeFunctionIdx++)
        {
            const uint32_t basisFunctionIdx = basisFunctionIndexer.getBasisFunctionIndex(elementIdx, shapeFunctionIdx);
            const Polynomial2D shapeFn = basisFunctionIndexer.getShapeFunction(elementIdx, shapeFunctionIdx);
            res += coefficients(basisFunctionIdx) * detA * integrateOverReferenceElement(shapeFn, element.getElementType());
        }
    }
    return res;
}
} // namespace fem
