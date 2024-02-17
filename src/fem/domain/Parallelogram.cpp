#include "fem/domain/Parallelogram.hpp"

namespace fem
{
Parallelogram::Parallelogram(const Node& n1, const Node& n2, const Node& n3, const Node& n4)
    : m_nodes({n1, n2, n3, n4})
{
    assert(areParallel(getSide(0), getSide(2)));
    assert(areParallel(getSide(1), getSide(3)));
    assert(getReferenceElementMap().A.determinant() > 0);
}

AffineMap Parallelogram::getReferenceElementMap() const
{
    Matrix2mpq A;
    A.col(0) = m_nodes[1] - m_nodes[0];
    A.col(1) = m_nodes[3] - m_nodes[0];
    const AffineMap F(A, m_nodes[0]);

    Matrix2mpq B;
    B.col(0) = Vector2mpq(mpq_class("1/2"), 0);
    B.col(1) = Vector2mpq(0, mpq_class("1/2"));
    const Vector2mpq c(mpq_class("1/2"), mpq_class("1/2"));
    const AffineMap G(B, c);

    return compose(F, G);
}
} // namespace fem
