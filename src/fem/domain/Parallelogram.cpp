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

std::vector<std::unique_ptr<Element>> Parallelogram::subdivide(const Vector2mpq& x) const
{
    std::vector<std::unique_ptr<Element>> res{};
    const AffineMap F = getReferenceElementMap();
    const AffineMap Finv = F.inverse();
    const Vector2mpq xloc = Finv(x);
    const bool pointIsOutside = !isPointInsideElement(x, *this);
    const bool pointIsInInterior = abs(xloc(0)) != 1 && abs(xloc(1)) != 1;
    const bool pointIsOnSide = (abs(xloc(0)) == 1 && abs(xloc(1)) != 1) || (abs(xloc(1)) == 1 && abs(xloc(0)) != 1);
    if (pointIsOutside)
    {
        res.push_back(std::make_unique<Parallelogram>(*this));
    }
    else if (pointIsInInterior)
    {
        const mpq_class tx = (xloc(0) + 1) / 2;
        const mpq_class ty = (xloc(1) + 1) / 2;
        const Vector2mpq dx = m_nodes[1] - m_nodes[0];
        const Vector2mpq dy = m_nodes[3] - m_nodes[0];
        const Vector2mpq n0 = m_nodes[0] + tx * dx;
        const Vector2mpq n1 = m_nodes[1] + ty * dy;
        const Vector2mpq n2 = m_nodes[3] + tx * dx;
        const Vector2mpq n3 = m_nodes[0] + ty * dy;
        res.push_back(std::make_unique<Parallelogram>(n0, x, n3, m_nodes[0]));
        res.push_back(std::make_unique<Parallelogram>(n1, x, n0, m_nodes[1]));
        res.push_back(std::make_unique<Parallelogram>(n2, x, n1, m_nodes[2]));
        res.push_back(std::make_unique<Parallelogram>(n3, x, n2, m_nodes[3]));
    }
    else if (pointIsOnSide)
    {
        uint32_t sideIdx;
        if (xloc(1) == -1)
        {
            sideIdx = 0;
        }
        else if (xloc(0) == 1)
        {
            sideIdx = 1;
        }
        else if (xloc(1) == 1)
        {
            sideIdx = 2;
        }
        else
        {
            sideIdx = 3;
        }
        const Vector2mpq d = m_nodes[(sideIdx + 2) % 4] - m_nodes[(sideIdx + 1) % 4];
        const Vector2mpq n0 = x + d;
        res.push_back(std::make_unique<Parallelogram>(m_nodes[sideIdx], x, n0, m_nodes[(sideIdx + 3) % 4]));
        res.push_back(std::make_unique<Parallelogram>(n0, x, m_nodes[(sideIdx + 1) % 4], m_nodes[(sideIdx + 2) % 4]));
    }
    else
    {
        res.push_back(std::make_unique<Parallelogram>(*this));
    }
    return res;
}
} // namespace fem
