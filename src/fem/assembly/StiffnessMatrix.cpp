#include "fem/assembly/StiffnessMatrix.hpp"

namespace fem
{
MatrixXmpq assembleStiffnessMatrix(const BasisFunctionIndexer& basisFunctionIndexer, const ShapeFunctionFactory& shapeFunctionFactory)
{
    return assembleStiffnessMatrix(basisFunctionIndexer);
}

MatrixXmpq assembleStiffnessMatrix(const BasisFunctionIndexer& basisFunctionIndexer)
{
    const uint32_t numOfBasisFunctions = basisFunctionIndexer.getNumOfBasisFunctions();
    MatrixXmpq res(numOfBasisFunctions, numOfBasisFunctions);
    const Mesh& mesh = basisFunctionIndexer.getMesh();
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
                const Polynomial2D shapeFunction1Dx = basisFunctionIndexer.getShapeFunctionDerivative(elementIdx, shapeFunctionIdx1, 'x');
                const Polynomial2D shapeFunction1Dy = basisFunctionIndexer.getShapeFunctionDerivative(elementIdx, shapeFunctionIdx1, 'y');
                const Polynomial2D shapeFunction2Dx = basisFunctionIndexer.getShapeFunctionDerivative(elementIdx, shapeFunctionIdx2, 'x');
                const Polynomial2D shapeFunction2Dy = basisFunctionIndexer.getShapeFunctionDerivative(elementIdx, shapeFunctionIdx2, 'y');
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

MatrixXmpq extractSubStiffnessMatrix(const MatrixXmpq& stiffnessMatrix, const BasisFunctionIndexer& superBasisFunctionIndexer, const BasisFunctionIndexer& subBasisFunctionIndexer)
{
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
