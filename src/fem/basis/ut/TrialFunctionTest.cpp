#include <gtest/gtest.h>

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
    const Mesh mesh{nodes, elements};
    const ShapeFunctionFactory shapeFunctionFactory{};
    const BasisFunctionIndexer basisFunctionIndexer{mesh, 3, PolynomialSpaceType_Product};
    VectorXmpq coefficients;

    void SetUp()
    {
        coefficients = VectorXmpq(basisFunctionIndexer.getNumOfBasisFunctions());
        coefficients(basisFunctionIndexer.getBasisFunctionIndex(0, 2)) = mpq_class("-1/4");
        coefficients(basisFunctionIndexer.getBasisFunctionIndex(0, 4 + 4 + 1)) = 3;
        coefficients(basisFunctionIndexer.getBasisFunctionIndex(0, 4 + 4*2)) = 1;
        coefficients(basisFunctionIndexer.getBasisFunctionIndex(1, 2)) = mpq_class("1/2");
    }
};

TEST_F(TrialFunctionTest, EvaluateTrialFunction)
{
    EXPECT_EQ(evaluateTrialFunction(coefficients, basisFunctionIndexer, {mpq_class("1/2"), mpq_class("-1/5")}), mpq_class("-141/200"));
    EXPECT_EQ(evaluateTrialFunction(coefficients, basisFunctionIndexer, {mpq_class("-3/8"), mpq_class("1")}), mpq_class("1445/512"));
    EXPECT_EQ(evaluateTrialFunction(coefficients, basisFunctionIndexer, {mpq_class("0"), mpq_class("2")}), mpq_class("1/2"));

    EXPECT_EQ(evaluateTrialFunction(coefficients, basisFunctionIndexer, {mpq_class("1/2"), mpq_class("-1/5")}, shapeFunctionFactory), mpq_class("-141/200"));
    EXPECT_EQ(evaluateTrialFunction(coefficients, basisFunctionIndexer, {mpq_class("-3/8"), mpq_class("1")}, shapeFunctionFactory), mpq_class("1445/512"));
    EXPECT_EQ(evaluateTrialFunction(coefficients, basisFunctionIndexer, {mpq_class("0"), mpq_class("2")}, shapeFunctionFactory), mpq_class("1/2"));
}

TEST_F(TrialFunctionTest, IntegrateTrialFunction)
{
    EXPECT_EQ(integrateTrialFunction(coefficients, basisFunctionIndexer), mpq_class("55/36") + mpq_class("1/12"));
    EXPECT_EQ(integrateTrialFunction(coefficients, basisFunctionIndexer, shapeFunctionFactory), mpq_class("55/36") + mpq_class("1/12"));

    coefficients = VectorXmpq(basisFunctionIndexer.getNumOfBasisFunctions());
    coefficients(0) = 1;
    coefficients(1) = 1;
    coefficients(2) = 1;
    coefficients(3) = 1;
    coefficients(4) = 1;
    EXPECT_EQ(integrateTrialFunction(coefficients, basisFunctionIndexer), mpq_class("5"));
    EXPECT_EQ(integrateTrialFunction(coefficients, basisFunctionIndexer, shapeFunctionFactory), mpq_class("5"));
}

TEST_F(TrialFunctionTest, NormalizeTrialFunction)
{
    normalizeTrialFunction(coefficients, basisFunctionIndexer, shapeFunctionFactory);
    EXPECT_EQ(integrateTrialFunction(coefficients, basisFunctionIndexer), 0);
    EXPECT_EQ(integrateTrialFunction(coefficients, basisFunctionIndexer, shapeFunctionFactory), 0);
}
} // namespace fem::ut
