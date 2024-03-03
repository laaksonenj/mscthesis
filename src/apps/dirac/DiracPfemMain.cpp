#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

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
    Timer timer;
    timer.start();

    Arguments args(argc, argv);

    const std::string meshFilename = args.getValue<std::string>("mesh-file");
    const uint32_t p_max = args.getValue<int>("p");
    const PolynomialSpaceType polynomialSpaceType = args.getValue<PolynomialSpaceType>("polynomial-space");
    const Vector2mpq x_0 = args.getValue<Vector2mpq>("dirac-point");
    const fs::path outputFilepath = fs::path(args.getValue<std::string>("output-file"));
    const uint32_t precision = args.getValue<uint32_t>("precision");

    std::cout << "Arguments:" << std::endl;
    std::cout << "--mesh-file " << meshFilename << std::endl;
    std::cout << "--p " << p_max << std::endl;
    std::cout << "--polynomial-space " << polynomialSpaceType << std::endl;
    std::cout << "--dirac-point " << x_0(0) << " " << x_0(1) << std::endl;
    std::cout << "--output-file " << outputFilepath << std::endl;
    std::cout << "--precision " << precision << std::endl;
    std::cout << std::endl;

    mpf_set_default_prec(precision);
    std::cout << "Precision of mpf_class: " << mpf_get_default_prec() << std::endl;
    std::cout << "Number of OpenMP threads: " << omp_get_max_threads() << std::endl;
    std::cout << std::endl;

    fs::create_directories(outputFilepath.parent_path());
    std::ofstream outputFile(outputFilepath);
    if (!outputFile.is_open())
    {
        std::cout << "Failed to open file " << outputFilepath << std::endl;
        return 1;
    }
    outputFile << std::setprecision(16);

    const auto mesh = std::make_shared<Mesh>(createMeshFromFile(meshFilename));
    const FemContext ctx(mesh, p_max, polynomialSpaceType);

    ShapeFunctionFactory shapeFunctionFactory;
    if (mesh->containsQuadrilateral())
    {
        std::cout << "Precomputing quadrilateral shape functions... ";
        timer.lap();
        shapeFunctionFactory.createShapeFunctions(ElementType_Parallelogram, p_max);
        std::cout << "done in " << minutesSecondsMilliseconds(timer.lap()) << std::endl;
    }
    if (mesh->containsTriangle())
    {
        std::cout << "Precomputing triangle shape functions... ";
        timer.lap();
        shapeFunctionFactory.createShapeFunctions(ElementType_Triangle, p_max);
        std::cout << "done in " << minutesSecondsMilliseconds(timer.lap()) << std::endl;
    }

    ShapeFunctionEvaluator shapeFunctionEvaluator(shapeFunctionFactory);
    if (mesh->containsQuadrilateral())
    {
        std::cout << "Pre-evaluating quadrilateral shape functions... ";
        timer.lap();
        shapeFunctionEvaluator.preEvaluate(ElementType_Parallelogram, defaultGLTableQuad.getAbscissas());
        std::cout << "done in " << minutesSecondsMilliseconds(timer.lap()) << std::endl;
    }
    if (mesh->containsTriangle())
    {
        std::cout << "Pre-evaluating triangle shape functions... ";
        timer.lap();
        shapeFunctionEvaluator.preEvaluate(ElementType_Triangle, defaultGLTableTri.getAbscissas());
        std::cout << "done in " << minutesSecondsMilliseconds(timer.lap()) << std::endl;
    }

    const auto exact = getNormalizedGreensFunction(x_0, *mesh);
    const auto grad_exact = getGreensFunctionGradient(x_0);

    std::cout << "Assembling stiffness matrix... ";
    timer.lap();
    const MatrixXmpq stiffnessMatrix = assembleStiffnessMatrix(ctx, shapeFunctionFactory);
    std::cout << "done in " << minutesSecondsMilliseconds(timer.lap()) << std::endl;

    std::cout << "Assembling Dirac load vector... ";
    timer.lap();
    const VectorXmpq diracLoadVector = assembleDiracLoadVector(ctx, x_0, shapeFunctionFactory);
    std::cout << "done in " << minutesSecondsMilliseconds(timer.lap()) << std::endl;

    std::cout << "Assembling Neumann load vector... ";
    timer.lap();
    const VectorXmpq neumannLoadVector = assembleNeumannLoadVector(ctx, grad_exact, shapeFunctionFactory);
    std::cout << "done in " << minutesSecondsMilliseconds(timer.lap()) << std::endl;

    const VectorXmpq loadVector = diracLoadVector + neumannLoadVector;

    std::cout << std::endl;
    std::cout << "--------------------------------------------------------------" << std::endl;
    for (int p = 1; p <= p_max; p++)
    {
        std::cout << "p=" << p << ":" << std::endl;

        const FemContext subCtx(ctx.mesh, p, ctx.polynomialSpaceType);

        std::cout << "Extracting system of equations... ";
        timer.lap();
        const MatrixXmpq subStiffnessMatrix = extractSubStiffnessMatrix(ctx, stiffnessMatrix, p);
        const VectorXmpq subLoadVector = extractSubLoadVector(ctx, loadVector, p);
        const uint32_t dim = subLoadVector.size();
        const MatrixXmpq A_mpq = subStiffnessMatrix.block(1, 1, dim-1, dim-1);
        const VectorXmpq b_mpq = subLoadVector.segment(1, dim-1);
        const MatrixXmpf A_mpf = A_mpq.cast<mpf_class>();
        const VectorXmpf b_mpf = b_mpq.cast<mpf_class>();
        std::cout << "done in " << minutesSecondsMilliseconds(timer.lap()) << std::endl;

        std::cout << "Solving system of equations... ";
        timer.lap();
        const VectorXmpf x_mpf = A_mpf.llt().solve(b_mpf);
        const VectorXmpq x_mpq = x_mpf.cast<mpq_class>();
        VectorXmpq coeffs(dim);
        coeffs(0) = 0;
        coeffs.segment(1, dim-1) = x_mpq;
        std::cout << "done in " << minutesSecondsMilliseconds(timer.lap()) << std::endl;

        std::cout << "Normalizing solution... ";
        timer.lap();
        normalizeTrialFunction(subCtx, coeffs, shapeFunctionFactory);
        std::cout << "done in " << minutesSecondsMilliseconds(timer.lap()) << std::endl;

        mpq_class squaredL2error = 0;
        for (int elementIdx = 0; elementIdx < mesh->getNumOfElements(); elementIdx++)
        {
            std::cout << "Computing L2 error over element " << elementIdx << "... ";
            timer.lap();
            const mpq_class err = computeSquaredL2ErrorOverElement(subCtx, coeffs, exact, elementIdx, x_0, shapeFunctionEvaluator);
            std::cout << "done in " << minutesSecondsMilliseconds(timer.lap()) << std::endl;
            squaredL2error += err;
        }
        const mpf_class L2error = sqrt(mpf_class(squaredL2error));
        outputFile << p << " " << L2error << std::endl;
        std::cout << "L2 error: " << L2error << std::endl;
        std::cout << "--------------------------------------------------------------" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Wall time: " << minutesSecondsMilliseconds(timer.getElapsedTime()) << std::endl;

    return 0;
}
