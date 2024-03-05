#pragma once

#include <map>
#include <string>

#include "fem/multiprecision/Types.hpp"

namespace fem
{
class LinearSolver
{
public:
    enum Method
    {
        PartialPivLU,
        ColPivHouseholderQR,
        LLT,
        LDLT,
        BDCSVD
    };

public:
    explicit LinearSolver(Method method);

    VectorXmpq solve(const MatrixXmpq& A, const VectorXmpq& b);
    mpf_class getRelativeError() const { return m_relativeError; }

private:
    Method m_method;
    mpf_class m_relativeError;
};

inline const std::map<LinearSolver::Method, std::string> linearSolverMethodCliNames{
    {LinearSolver::PartialPivLU, "partial-piv-lu"},
    {LinearSolver::ColPivHouseholderQR, "col-piv-householder-qr"},
    {LinearSolver::LLT, "llt"},
    {LinearSolver::LDLT, "ldlt"},
    {LinearSolver::BDCSVD, "bdcsvd"}
};
} // namespace fem
