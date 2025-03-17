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
    int gridpoints;               // number of gridpoints
    int niter;                    // number of iterations for the solver
    bool BS_solver_verbosity;     // verbosity for BS solver
    bool BS_enable_matching;      // whether to enable macthing to asymptotic
                                  // expressions
    double PSC;                   // central value of the conformal factor
    double OMC;                   // central value of the lapse
    double central_amplitude_CSF; // central amplitude
    double scalar_mass;           // mu (usually = 1)
    double phi4_coeff;            // quartic term in the potential
    bool solitonic;         // whether potential should be of solitonic kind
    double sigma_solitonic; // self-interaction term for the solitonic potential
    double phase;           // phase offset
    bool antiboson; // whether one star should have an oppositive frequency to
                    // its companion in a binary
    double BS_separation;       // separation along the x-axis
    double BS_impact_parameter; // separation along the y-axis
    double BS_rapidity;         // rapidity (atanh(v))
    double mass;                // mass of a BS
    double mass_ratio;          // mass ratio
    double radius_width1;       // R_A for unequal-mass fix
    double radius_width2;       // R_B for unequal-mass fix
    int conformal_factor_power; // power for the conformal factor expression in
                                // the unequal-mass fix
    std::array<double, CH_SPACEDIM>
        star_centre; // coordinates of the centre of the star
};

#endif /* BOSONSTARPARAMS_HPP_ */
