/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLEPRECONDITIONER_HPP
#define DIAGNOSTICVARIABLEPRECONDITIONER_HPP

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.H"
#include "MatterCCZ4.hpp"
#include "simd.hpp"

/*
 * Sets a diagnostic field variable to zero outside a set radius (a_radius)
 */

template <class matter_t> class DiagnosticVariablePreconditioner
{
  public:
    template <class data_t>
    using Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

    //! Constructor
    DiagnosticVariablePreconditioner(
        const matter_t &a_matter, const double dx, const double a_L,
        const std::array<double, CH_SPACEDIM> a_centre, const int a_c_var1 = -1,
        const int a_c_var2 = -1, const double a_radius = 0.)
        : m_matter(a_matter), m_deriv(dx), m_dx(dx), m_L(a_L),
          m_centre(a_centre), m_c_var1(a_c_var1), m_c_var2(a_c_var2),
          m_radius(a_radius)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx, m_centre);
        const auto vars = current_cell.template load_vars<Vars>();

        data_t r = coords.get_radius();
        data_t r_max = m_radius;
        data_t heaviside = 0.5 * (1. - abs(r - r_max) / (r - r_max));

        auto original_val1 = current_cell.load_vars(m_c_var1);
        data_t output_val1 = original_val1 * heaviside;
        current_cell.store_vars(output_val1, m_c_var1);

        auto original_val2 = current_cell.load_vars(m_c_var2);
        data_t output_val2 = original_val2 * heaviside;
        current_cell.store_vars(output_val2, m_c_var2);
    }

  protected:
    const double m_L, m_dx, m_radius;
    const matter_t &m_matter;
    FourthOrderDerivatives m_deriv;
    const int m_c_var1, m_c_var2; // var enum
    const std::array<double, CH_SPACEDIM> m_centre;
};

#endif /* DIAGNOSTICVARIABLEPRECONDITIONER_HPP */
