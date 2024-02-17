#include "fem/domain/Triangle.hpp"

namespace fem
{
Triangle::Triangle(const Node& n1, const Node& n2, const Node& n3)
    : m_nodes({n1, n2, n3})
{
    assert(getReferenceElementMap().A.determinant() > 0);
}

AffineMap Triangle::getReferenceElementMap() const
{
    Matrix2mpq A;
    A.col(0) = m_nodes[1] - m_nodes[0];
    A.col(1) = m_nodes[2] - m_nodes[0];
    return AffineMap(A, m_nodes[0]);
}
} // namespace fem
