#include "apps/common/LinearSolver.hpp"

#include <cassert>

#include <Eigen/Dense>

namespace fem
{
namespace
{
mpf_class computeRelativeError(const MatrixXmpq& A, const VectorXmpq& x, const VectorXmpq& b)
{
    const VectorXmpq err = A*x - b;
    return sqrt(mpf_class(err.squaredNorm())) / sqrt(mpf_class(b.squaredNorm()));
}
} // namespace

LinearSolver::LinearSolver(Method method)
    : m_method(method)
    , m_relativeError(-1)
{
}

VectorXmpq LinearSolver::solve(const MatrixXmpq& A, const VectorXmpq& b)
{
    VectorXmpq res;
    const Eigen::MatrixXd A_d = A.unaryExpr([](const mpq_class& elem) -> double { return elem.get_d(); });
    const Eigen::VectorXd b_d = b.unaryExpr([](const mpq_class& elem) -> double { return elem.get_d(); });
    if (m_method == PartialPivLU)
    {
        res = A.partialPivLu().solve(b); // becomes extremely slow very quickly so use mainly for validation etc.
    }
    else if (m_method == ColPivHouseholderQR)
    {
        const Eigen::VectorXd x_d = A_d.colPivHouseholderQr().solve(b_d);
        res = x_d.cast<mpq_class>();
    }
    else if (m_method == LLT)
    {
        const Eigen::VectorXd x_d = A_d.llt().solve(b_d);
        res = x_d.cast<mpq_class>();
    }
    else if (m_method == LDLT)
    {
        const Eigen::VectorXd x_d = A_d.ldlt().solve(b_d);
        res = x_d.cast<mpq_class>();
    }
    else if (m_method == BDCSVD)
    {
        const Eigen::VectorXd x_d = A_d.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_d);
        res = x_d.cast<mpq_class>();
    }
    else
    {
        assert(false && "Unknown method");
    }
    m_relativeError = computeRelativeError(A, res, b);
    return res;
}
} // namespace fem
