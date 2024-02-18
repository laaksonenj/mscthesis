#include <cstdint>
#include <iostream>
#include <sstream>

#include "apps/common/Arguments.hpp"
#include "apps/refdata/Utils.hpp"
#include "fem/assembly/StiffnessMatrix.hpp"
#include "fem/basis/PolynomialSpaceType.hpp"
#include "fem/domain/Mesh.hpp"

using namespace fem;

int main(int argc, char* argv[])
{
    Arguments args(argc, argv);

    const std::string meshFilename = args.getValue<std::string>("mesh-file");
    const uint32_t p = args.getValue<int>("p");
    const PolynomialSpaceType polynomialSpaceType = args.getValue<PolynomialSpaceType>("polynomial-space");
    const std::string varName = args.getValue<std::string>("variable-name");
    
    const Mesh mesh = createMeshFromFile(meshFilename);
    const BasisFunctionIndexer basisFunctionIndexer(mesh, p, polynomialSpaceType);
    const MatrixXmpq stiffnessMatrix = assembleStiffnessMatrix(basisFunctionIndexer);

    std::stringstream ss;
    ss << generateCommonHeader();
    ss << namespaceBegin();
    ss << generateVariableDefinition(stiffnessMatrix, varName);
    ss << namespaceEnd();
    std::cout << ss.str();

    return 0;
}
