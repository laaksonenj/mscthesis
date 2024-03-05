#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <gmpxx.h>
#include <omp.h>

#include "apps/common/Arguments.hpp"
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
    const fs::path outputFilepath = fs::path(args.getValue<std::string>("output-file"));

    std::cout << "Arguments:" << std::endl;
    std::cout << "--mesh-file " << meshFilename << std::endl;
    std::cout << "--p " << p_max << std::endl;
    std::cout << "--polynomial-space " << polynomialSpaceType << std::endl;
    std::cout << "--dirac-point " << x_0(0) << " " << x_0(1) << std::endl;
    std::cout << "--output-file " << outputFilepath << std::endl;
    std::cout << std::endl;

    std::cout << "Number of OpenMP threads: " << omp_get_max_threads() << std::endl;
    std::cout << std::endl;

    std::ofstream outputFile;
    if (!outputFilepath.empty())
    {
        if (!outputFilepath.parent_path().empty())
        {
            fs::create_directories(outputFilepath.parent_path());
        }
        outputFile.open(outputFilepath);
        if (!outputFile.is_open())
        {
            std::cout << "Failed to create output file " << outputFilepath << std::endl;
            return 1;
        }
        outputFile << std::setprecision(16);
    }

    const auto mesh = std::make_shared<Mesh>(createMeshFromFile(meshFilename));
    const FemContext ctx(mesh, p_max, polynomialSpaceType);

    Timer timer;

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
        shapeFunctionEvaluator.preEvaluate(ElementType_Parallelogram, defaultGLTableQuad.getAbscissas());
        timer.stop();
    }
    if (mesh->containsTriangle())
    {
        timer.start("Pre-evaluating triangle shape functions... ");
        shapeFunctionEvaluator.preEvaluate(ElementType_Triangle, defaultGLTableTri.getAbscissas());
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
        const MatrixXmpq A_mpq = subStiffnessMatrix.block(1, 1, dim-1, dim-1);
        const VectorXmpq b_mpq = subLoadVector.segment(1, dim-1);
        const Eigen::MatrixXd A_d = A_mpq.unaryExpr([](const mpq_class& elem) -> double { return elem.get_d(); });
        const Eigen::VectorXd b_d = b_mpq.unaryExpr([](const mpq_class& elem) -> double { return elem.get_d(); });
        timer.stop();

        timer.start("Solving system of equations... ");
        const Eigen::VectorXd x_d = A_d.llt().solve(b_d);
        const VectorXmpq x_mpq = x_d.cast<mpq_class>();
        VectorXmpq coeffs(dim);
        coeffs(0) = 0;
        coeffs.segment(1, dim-1) = x_mpq;
        timer.stop();

        timer.start("Normalizing solution... ");
        normalizeTrialFunction(subCtx, coeffs, shapeFunctionFactory);
        timer.stop();

        mpq_class squaredL2error = 0;
        for (int elementIdx = 0; elementIdx < mesh->getNumOfElements(); elementIdx++)
        {
            timer.start("Computing L2 error over element " + std::to_string(elementIdx) + "... ");
            const mpq_class err = computeSquaredL2ErrorOverElement(subCtx, coeffs, exact, elementIdx, x_0, shapeFunctionEvaluator);
            timer.stop();
            squaredL2error += err;
        }
        const mpf_class L2error = sqrt(mpf_class(squaredL2error));
        outputFile << p << " " << L2error << std::endl;
        std::cout << "L2 error: " << L2error << std::endl;
        std::cout << "--------------------------------------------------------------" << std::endl;
    }

    return 0;
}
