#pragma once

#include <Eigen/Dense>
#include <gmpxx.h>

namespace fem
{
using Vector2mpq = Eigen::Matrix<mpq_class, 2, 1>;
using Vector3mpq = Eigen::Matrix<mpq_class, 3, 1>;
using Matrix2mpq = Eigen::Matrix<mpq_class, 2, 2>;
using Matrix3mpq = Eigen::Matrix<mpq_class, 3, 3>;
using MatrixXmpq = Eigen::Matrix<mpq_class, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXmpq = Eigen::Matrix<mpq_class, Eigen::Dynamic, 1>;

using Vector2mpf = Eigen::Matrix<mpf_class, 2, 1>;
using Vector3mpf = Eigen::Matrix<mpf_class, 3, 1>;
using Matrix2mpf = Eigen::Matrix<mpf_class, 2, 2>;
using Matrix3mpf = Eigen::Matrix<mpf_class, 3, 3>;
using MatrixXmpf = Eigen::Matrix<mpf_class, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXmpf = Eigen::Matrix<mpf_class, Eigen::Dynamic, 1>;
} // namespace fem
