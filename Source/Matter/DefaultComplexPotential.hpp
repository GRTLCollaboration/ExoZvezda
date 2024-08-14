/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DEFAULTCOMPLEXPOTENTIAL_HPP_
#define DEFAULTCOMPLEXPOTENTIAL_HPP_

#include "Tensor.hpp"
#include "simd.hpp"

class DefaultComplexPotential
{
  public:
    //! The constructor
    DefaultComplexPotential() {}

    //! Set the potential function for the complex scalar field here to zero
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_modulus_phi_squared,
                           data_t &dVdmodulus_phi_squared,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at |phi|^2
        V_of_modulus_phi_squared = 0.0;

        // The potential gradient at |phi|^2
        dVdmodulus_phi_squared = 0.0;
    }
};

#endif /* DEFAULTCOMPLEXPOTENTIAL_HPP_ */
