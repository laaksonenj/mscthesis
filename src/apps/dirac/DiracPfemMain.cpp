#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <omp.h>

#include "apps/common/Arguments.hpp"
#include "apps/common/LinearSolver.hpp"
#include "apps/common/Timer.hpp"
#include "apps/dirac/utils/GreensFunction.hpp"
#include "apps/dirac/utils/L2Error.hpp"
#include "fem/basis/ShapeFunctionEvaluator.hpp"
#include "fem/basis/ShapeFunctionFactory.hpp"
#include "fem/basis/TrialFunction.hpp"
#include "fem/assembly/LoadVector.hpp"
#include "fem/assembly/StiffnessMatrix.hpp"
#include "fem/math/Quadrature.hpp"

using namespace fem;
namespace fs = std::filesystem;

int main(int argc, char* argv[])
{
    Arguments args(argc, argv);

    const std::string meshFilename = args.getValue<std::string>("mesh-file");
    const uint32_t p_max = args.getValue<int>("p");
    const PolynomialSpaceType polynomialSpaceType = args.getValue<PolynomialSpaceType>("polynomial-space");
    const Vector2mpq x_0 = args.getValue<Vector2mpq>("dirac-point");
    const fs::path outputDirpath = fs::path(args.getValue<std::string>("output-dir"));
    const LinearSolver::Method linearSolverMethod = args.getValue<LinearSolver::Method>("linear-solver");

    std::cout << "Arguments:" << std::endl;
    std::cout << "--mesh-file " << meshFilename << std::endl;
    std::cout << "--p " << p_max << std::endl;
    std::cout << "--polynomial-space " << polynomialSpaceType << std::endl;
    std::cout << "--dirac-point " << x_0(0) << " " << x_0(1) << std::endl;
    std::cout << "--output-dir " << outputDirpath << std::endl;
    std::cout << "--linear-solver " << linearSolverMethodCliNames.at(linearSolverMethod) << std::endl;
    std::cout << std::endl;

    std::cout << "Number of OpenMP threads: " << omp_get_max_threads() << std::endl;
    std::cout << std::endl;

    Timer timer;
    const auto mesh = std::make_shared<Mesh>(createMeshFromFile(meshFilename));
    const FemContext ctx(mesh, p_max, polynomialSpaceType);
    LinearSolver linearSolver(linearSolverMethod);

    std::ofstream globalErrorOutputFile;
    std::vector<std::ofstream> elementErrorOutputFiles(mesh->getNumOfElements());
    if (!outputDirpath.empty())
    {
        fs::create_directories(outputDirpath);
        const fs::path globalErrorFilepath = outputDirpath / "L2error.txt";
        globalErrorOutputFile.open(globalErrorFilepath);
        if (!globalErrorOutputFile.is_open())
        {
            std::cout << "Failed to create output file " << globalErrorFilepath  << std::endl;
            return 1;
        }
        globalErrorOutputFile << std::setprecision(16);
        for (int elementIdx = 0; elementIdx < mesh->getNumOfElements(); elementIdx++)
        {
            auto& ofs = elementErrorOutputFiles.at(elementIdx);
            const std::string elementErrorFilename = "L2error_element" + std::to_string(elementIdx) + ".txt";
            const fs::path elementErrorFilepath = outputDirpath / elementErrorFilename;
            ofs.open(elementErrorFilepath);
            if (!ofs.is_open())
            {
                std::cout << "Failed to create output file " << elementErrorFilepath  << std::endl;
                return 1;
            }
            ofs << std::setprecision(16);
        }
    }

    ShapeFunctionFactory shapeFunctionFactory;
    if (mesh->containsQuadrilateral())
    {
        timer.start("Precomputing quadrilateral shape functions... ");
        shapeFunctionFactory.createShapeFunctions(ElementType_Parallelogram, p_max);
        timer.stop();
    }
    if (mesh->containsTriangle())
    {
        timer.start("Precomputing triangle shape functions... ");
        shapeFunctionFactory.createShapeFunctions(ElementType_Triangle, p_max);
        timer.stop();
    }

    ShapeFunctionEvaluator shapeFunctionEvaluator(shapeFunctionFactory);
    if (mesh->containsQuadrilateral())
    {
        timer.start("Pre-evaluating quadrilateral shape functions... ");
        const auto preEvaluationPoints = getShapeFunctionEvaluationPointsForL2Error(ElementType_Parallelogram, *mesh, x_0);
        shapeFunctionEvaluator.preEvaluate(ElementType_Parallelogram, preEvaluationPoints);
        timer.stop();
    }
    if (mesh->containsTriangle())
    {
        timer.start("Pre-evaluating triangle shape functions... ");
        const auto preEvaluationPoints = getShapeFunctionEvaluationPointsForL2Error(ElementType_Triangle, *mesh, x_0);
        shapeFunctionEvaluator.preEvaluate(ElementType_Triangle, preEvaluationPoints);
        timer.stop();
    }

    const auto exact = getNormalizedGreensFunction(x_0, *mesh);
    const auto grad_exact = getGreensFunctionGradient(x_0);

    timer.start("Assembling stiffness matrix... ");
    const MatrixXmpq stiffnessMatrix = assembleStiffnessMatrix(ctx, shapeFunctionFactory);
    timer.stop();

    timer.start("Assembling Dirac load vector... ");
    const VectorXmpq diracLoadVector = assembleDiracLoadVector(ctx, x_0, shapeFunctionFactory);
    timer.stop();

    timer.start("Assembling Neumann load vector... ");
    const VectorXmpq neumannLoadVector = assembleNeumannLoadVector(ctx, grad_exact, shapeFunctionFactory);
    timer.stop();

    const VectorXmpq loadVector = diracLoadVector + neumannLoadVector;

    std::cout << std::endl;
    std::cout << "--------------------------------------------------------------" << std::endl;
    for (int p = 1; p <= p_max; p++)
    {
        std::cout << "p=" << p << ":" << std::endl;

        const FemContext subCtx(ctx.mesh, p, ctx.polynomialSpaceType);

        timer.start("Extracting system of equations... ");
        const MatrixXmpq subStiffnessMatrix = extractSubStiffnessMatrix(ctx, stiffnessMatrix, p);
        const VectorXmpq subLoadVector = extractSubLoadVector(ctx, loadVector, p);
        const uint32_t dim = subLoadVector.size();
        const MatrixXmpq A = subStiffnessMatrix.block(1, 1, dim-1, dim-1);
        const VectorXmpq b = subLoadVector.segment(1, dim-1);
        timer.stop();

        timer.start("Solving system of equations... ");
        const VectorXmpq x = linearSolver.solve(A, b);
        timer.stop();

        std::cout << "Relative error of solution due to floating-point: " << linearSolver.getRelativeError() << std::endl;

        VectorXmpq coeffs(dim);
        coeffs(0) = 0;
        coeffs.segment(1, dim-1) = x;

        timer.start("Normalizing solution... ");
        normalizeTrialFunction(subCtx, coeffs, shapeFunctionFactory);
        timer.stop();

        mpq_class squaredL2error = 0;
        for (int elementIdx = 0; elementIdx < mesh->getNumOfElements(); elementIdx++)
        {
            timer.start("Computing L2 error over element " + std::to_string(elementIdx) + "... ");
            const mpq_class err = computeSquaredL2ErrorOverElement(subCtx, coeffs, exact, elementIdx, x_0, shapeFunctionEvaluator);
            timer.stop();
            elementErrorOutputFiles.at(elementIdx) << p << " " << sqrt(mpf_class(err)) << std::endl;
            squaredL2error += err;
        }
        const mpf_class L2error = sqrt(mpf_class(squaredL2error));
        globalErrorOutputFile << p << " " << L2error << std::endl;
        std::cout << "L2 error: " << L2error << std::endl;
        std::cout << "--------------------------------------------------------------" << std::endl;
    }

    return 0;
}
