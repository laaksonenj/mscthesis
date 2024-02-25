#include <cstdint>
#include <iostream>
#include <sstream>

#include "apps/common/Arguments.hpp"
#include "apps/refdata/Utils.hpp"
#include "fem/assembly/DiracLoadVector.hpp"
#include "fem/basis/PolynomialSpaceType.hpp"
#include "fem/domain/Mesh.hpp"

using namespace fem;

int main(int argc, char* argv[])
{
    Arguments args(argc, argv);

    const std::string meshFilename = args.getValue<std::string>("mesh-file");
    const uint32_t p = args.getValue<int>("p");
    const PolynomialSpaceType polynomialSpaceType = args.getValue<PolynomialSpaceType>("polynomial-space");
    const Vector2mpq x_0 = args.getValue<Vector2mpq>("dirac-point");
    const std::string varName = args.getValue<std::string>("variable-name");

    const FemContext ctx(std::make_shared<Mesh>(createMeshFromFile(meshFilename)), p, polynomialSpaceType);
    const VectorXmpq diracLoadVector = assembleDiracLoadVector(ctx, x_0);

    std::stringstream ss;
    ss << generateCommonHeader();
    ss << namespaceBegin();
    ss << generateVariableDefinition(diracLoadVector, varName);
    ss << namespaceEnd();
    std::cout << ss.str();

    return 0;
}
