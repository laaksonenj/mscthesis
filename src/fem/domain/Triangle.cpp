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

std::vector<std::unique_ptr<Element>> Triangle::subdivideImpl(const Vector2mpq& x) const
{
    std::vector<std::unique_ptr<Element>> res{};
    const AffineMap F = getReferenceElementMap();
    const AffineMap Finv = F.inverse();
    const Vector2mpq xloc = Finv(x);
    const bool pointIsInside = xloc(0) > 0 && xloc(1) > 0 && xloc(0) + xloc(1) < 1;
    const bool pointIsOnSide = (xloc(0) == 0 && xloc(1) > 0 && xloc(1) < 1)
                               || (xloc(1) == 0 && xloc(0) > 0 && xloc(0) < 1)
                               || (xloc(0) + xloc(1) == 1 && xloc(0) != 1 && xloc(1) != 1);
    if (pointIsInside)
    {
        res.push_back(std::make_unique<Triangle>(m_nodes[0], x, m_nodes[2]));
        res.push_back(std::make_unique<Triangle>(m_nodes[1], x, m_nodes[0]));
        res.push_back(std::make_unique<Triangle>(m_nodes[2], x, m_nodes[1]));
    }
    else if (pointIsOnSide)
    {
        uint32_t sideIdx;
        if (xloc(1) == 0)
        {
            sideIdx = 0;
        }
        else if (xloc(0) == 0)
        {
            sideIdx = 2;
        }
        else
        {
            sideIdx = 1;
        }
        res.push_back(std::make_unique<Triangle>(m_nodes[sideIdx], x, m_nodes[(sideIdx + 2) % 3]));
        res.push_back(std::make_unique<Triangle>(m_nodes[(sideIdx + 2) % 3], x, m_nodes[(sideIdx + 1) % 3]));
    }
    else
    {
        res.push_back(std::make_unique<Triangle>(*this));
    }
    return res;
}
} // namespace fem
