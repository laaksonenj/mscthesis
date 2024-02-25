#include <cstdint>
#include <iostream>
#include <sstream>

#include "apps/common/Arguments.hpp"
#include "apps/refdata/Utils.hpp"
#include "fem/assembly/NeumannLoadVector.hpp"
#include "fem/assembly/ut/ReferenceBoundaryFunctions.hpp"
#include "fem/basis/PolynomialSpaceType.hpp"
#include "fem/domain/Mesh.hpp"

using namespace fem;

int main(int argc, char* argv[])
{
    Arguments args(argc, argv);

    const std::string meshFilename = args.getValue<std::string>("mesh-file");
    const uint32_t p = args.getValue<int>("p");
    const PolynomialSpaceType polynomialSpaceType = args.getValue<PolynomialSpaceType>("polynomial-space");
    const uint32_t boundaryFunctionIdx = args.getValue<int>("boundary-function-idx");
    const uint32_t elementIdx = args.getValue<int>("element-idx");
    const uint32_t localSideIdx = args.getValue<int>("local-side-idx");
    const std::string varName = args.getValue<std::string>("variable-name");

    const FemContext ctx(std::make_shared<Mesh>(createMeshFromFile(meshFilename)), p, polynomialSpaceType);
    const VectorXmpq neumannLoadVector = assembleNeumannLoadVector(ctx, ut::referenceBoundaryFunctions.at(boundaryFunctionIdx), elementIdx, localSideIdx);

    std::stringstream ss;
    ss << generateCommonHeader();
    ss << namespaceBegin();
    ss << generateVariableDefinition(neumannLoadVector, varName);
    ss << namespaceEnd();
    std::cout << ss.str();

    return 0;
}
