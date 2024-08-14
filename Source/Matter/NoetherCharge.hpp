/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef NOETHERCHARGE_HPP_
#define NOETHERCHARGE_HPP_

#include "ADMConformalVars.hpp" // needed for CCz4 and matter variables
#include "Cell.hpp"
#include "ComplexScalarField.hpp"
#include "Coordinates.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"

//! Calculates the Noether Charge integrand values and the modulus of the
//! complex scalar field on the grid
class NoetherCharge
{
    // Need matter variables and chi
    template <class data_t>
    using ADMVars = ADMConformalVars::VarsNoGauge<data_t>;
    template <class data_t>
    using MatterVars = ComplexScalarField<>::Vars<data_t>;

  public:
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars locally
        const auto adm_vars = current_cell.template load_vars<ADMVars>();
        const auto matter_vars = current_cell.template load_vars<MatterVars>();

        // calculate Noether charge
        data_t N =
            pow(adm_vars.chi, -1.5) * (matter_vars.phi_Im * matter_vars.Pi_Re -
                                       matter_vars.phi_Re * matter_vars.Pi_Im);

        data_t mod_phi = sqrt(matter_vars.phi_Re * matter_vars.phi_Re +
                              matter_vars.phi_Im * matter_vars.phi_Im);

        current_cell.store_vars(N, c_N);
        current_cell.store_vars(mod_phi, c_mod_phi);
    }
};

#endif /* NOETHERCHARGE_HPP_ */
