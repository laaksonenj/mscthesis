#include <gtest/gtest.h>

#include "fem/basis/BasisFunctionIndexer.hpp"
#include "fem/basis/TrialFunction.hpp"

namespace fem::ut
{
class TrialFunctionTest : public testing::Test
{
protected:
    const std::vector<Node> nodes{
        {-1, -1},
        {1, -1},
        {1, 1},
        {-1, 1},
        {0, 2}
    };
    const std::vector<std::vector<Mesh::NodeIndex>> elements{
        {0, 1, 2, 3},
        {3, 2, 4}
    };
    const std::shared_ptr<Mesh> mesh{new Mesh(nodes, elements)};
    const uint32_t p{3};
    const FemContext ctx{mesh, p, PolynomialSpaceType_Product};
    const BasisFunctionIndexer basisFunctionIndexer{ctx};
    ShapeFunctionFactory shapeFunctionFactory{};
    const ShapeFunctionEvaluator shapeFunctionEvaluator{shapeFunctionFactory};
    VectorXmpq coefficients;

    void SetUp()
    {
        shapeFunctionFactory.createShapeFunctions(ElementType_Parallelogram, p);
        shapeFunctionFactory.createShapeFunctions(ElementType_Triangle, p);
        coefficients = VectorXmpq(basisFunctionIndexer.getNumOfBasisFunctions());
        coefficients(basisFunctionIndexer.getBasisFunctionIndex(0, 2)) = mpq_class("-1/4");
        coefficients(basisFunctionIndexer.getBasisFunctionIndex(0, 4 + 4 + 1)) = 3;
        coefficients(basisFunctionIndexer.getBasisFunctionIndex(0, 4 + 4*2)) = 1;
        coefficients(basisFunctionIndexer.getBasisFunctionIndex(1, 2)) = mpq_class("1/2");
    }
};

TEST_F(TrialFunctionTest, EvaluateTrialFunction)
{
    EXPECT_EQ(evaluateTrialFunction(ctx, coefficients, {mpq_class("1/2"), mpq_class("-1/5")}), mpq_class("-141/200"));
    EXPECT_EQ(evaluateTrialFunction(ctx, coefficients, {mpq_class("-3/8"), mpq_class("1")}), mpq_class("1445/512"));
    EXPECT_EQ(evaluateTrialFunction(ctx, coefficients, {mpq_class("0"), mpq_class("2")}), mpq_class("1/2"));

    EXPECT_EQ(evaluateTrialFunction(ctx, coefficients, {mpq_class("1/2"), mpq_class("-1/5")}, shapeFunctionEvaluator), mpq_class("-141/200"));
    EXPECT_EQ(evaluateTrialFunction(ctx, coefficients, {mpq_class("-3/8"), mpq_class("1")}, shapeFunctionEvaluator), mpq_class("1445/512"));
    EXPECT_EQ(evaluateTrialFunction(ctx, coefficients, {mpq_class("0"), mpq_class("2")}, shapeFunctionEvaluator), mpq_class("1/2"));
}

TEST_F(TrialFunctionTest, IntegrateTrialFunction)
{
    EXPECT_EQ(integrateTrialFunction(ctx, coefficients), mpq_class("55/36") + mpq_class("1/12"));
    EXPECT_EQ(integrateTrialFunction(ctx, coefficients, shapeFunctionFactory), mpq_class("55/36") + mpq_class("1/12"));

    coefficients = VectorXmpq(basisFunctionIndexer.getNumOfBasisFunctions());
    coefficients(0) = 1;
    coefficients(1) = 1;
    coefficients(2) = 1;
    coefficients(3) = 1;
    coefficients(4) = 1;
    EXPECT_EQ(integrateTrialFunction(ctx, coefficients), mpq_class("5"));
    EXPECT_EQ(integrateTrialFunction(ctx, coefficients, shapeFunctionFactory), mpq_class("5"));
}

TEST_F(TrialFunctionTest, NormalizeTrialFunction)
{
    normalizeTrialFunction(ctx, coefficients, shapeFunctionFactory);
    EXPECT_EQ(integrateTrialFunction(ctx, coefficients), 0);
    EXPECT_EQ(integrateTrialFunction(ctx, coefficients, shapeFunctionFactory), 0);
}
} // namespace fem::ut
