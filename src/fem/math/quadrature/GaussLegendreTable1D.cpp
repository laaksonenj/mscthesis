#include "fem/math/quadrature/GaussLegendreTable1D.hpp"

#include <cassert>

#include <gsl/gsl_integration.h>

namespace fem
{
GaussLegendreTable1D::GaussLegendreTable1D(uint32_t n)
{
    assert(n > 0);
    gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(n);
    double x, w;
    for (int i = 0; i < n; i++)
    {
        gsl_integration_glfixed_point(-1, 1, i, &x, &w, t);
        m_weights.push_back(mpq_class(w));
        m_abscissas.push_back(mpq_class(x));
    }
    gsl_integration_glfixed_table_free(t);
}
} // namespace fem
