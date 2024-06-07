/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>
#include <vector>

#ifndef ANGMOMFLUXPARAMS_HPP_
#define ANGMOMFLUXPARAMS_HPP_

//! A structure for the input params for the boson star
struct AngMomFlux_params_t
{
    int number_radii;
    bool do_flux_integration; // need to put this in specific level
    int extraction_level;
    int num_phi, num_theta;
    std::array<double, CH_SPACEDIM> centre;
    std::vector<double> radii;

    // need to sort out the field numbers... try mirens interval thing
};

#endif /* ANGMOMFLUXPARAMS_HPP_ */
