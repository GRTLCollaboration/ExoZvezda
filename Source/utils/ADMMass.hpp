/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMMASS_HPP_
#define ADMMASS_HPP_

#include "ADMConformalVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the ADM massa
class ADMMass
{
    // Use the variable definition in ADMConformalVars - only require the key
    // vars
    template <class data_t> using Vars = ADMConformalVars::VarsNoGauge<data_t>;

    template <class data_t>
    using Diff1Vars = ADMConformalVars::Diff2VarsNoGauge<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_centre;

  public:
    ADMMass(const std::array<double, CH_SPACEDIM> a_centre, double a_dx)
        : m_deriv(a_dx), m_dx(a_dx), m_centre(a_centre)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // Load the required vars
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Diff1Vars>(current_cell);

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);

        Coordinates<data_t> coords(current_cell, m_dx, m_centre);
        data_t r = coords.get_radius();

        Tensor<1, data_t> x = {coords.x, coords.y, coords.z};

        // Normal to the surface
        Tensor<1, data_t> dS = {coords.x / r, coords.y / r, coords.z / r};

        // ADM mass
        data_t Madm = 0.0;

        FOR4(i, j, k, l)
        {
            Madm += dS[i] / (16. * M_PI) * h_UU[j][k] * h_UU[i][l] *
                    (pow(vars.chi, -0.5) * (d1.h[l][k][j] - d1.h[j][k][l]) -
                     pow(vars.chi, -1.5) *
                         (vars.h[l][j] * d1.chi[j] - vars.h[j][k] * d1.chi[l]));
        }

        // Spin about z axis
        data_t Jadm = 0.0;

        // Note this is the levi civita symbol,
        const Tensor<3, double> epsilon = TensorAlgebra::epsilon();

        FOR3(i, j, k)
        {
            Jadm += -dS[i] / (8. * M_PI) * epsilon[2][j][k] *
                    pow(vars.chi, -1.5) * x[j] * vars.K *
                    TensorAlgebra::delta(i, k);

            FOR2(l, m)
            {
                Jadm += dS[i] / (8. * M_PI) * epsilon[2][j][k] * x[j] *
                        h_UU[i][l] * h_UU[k][m] * pow(vars.chi, -0.5) *
                        (vars.A[l][m] + 1.0 / 3.0 * vars.K * vars.h[l][m]);
            }
        }

        current_cell.store_vars(Madm, c_Madm);
        current_cell.store_vars(Jadm, c_Jadm);
    }
};

#endif /* ADMMASS_HPP_ */
