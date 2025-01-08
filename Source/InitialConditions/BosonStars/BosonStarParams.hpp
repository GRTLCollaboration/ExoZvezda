/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>

#ifndef BOSONSTARPARAMS_HPP_
#define BOSONSTARPARAMS_HPP_

//! A structure for the input params for the boson star
struct BosonStar_params_t
{
    int gridpoints;
    double PSC;
    double OMC;
    bool BS_solver_verbosity;
    bool BS_enable_matching;
    double central_amplitude_CSF;
    double scalar_mass;
    double phi4_coeff;
    bool solitonic;
    double sigma_solitonic; // 0.0499;
    double phase;
    double BS_separation;
    double mass;
    double BS_impact_parameter;
    double BS_rapidity;
    double mass_ratio;
    int n_power;
    int niter;
    double radius_width1;
    double radius_width2;
    int conformal_factor_power;
    bool antiboson;
    std::array<double, CH_SPACEDIM>
        star_centre; //!< coordinates of the centre of the star
};

#endif /* BOSONSTARPARAMS_HPP_ */
