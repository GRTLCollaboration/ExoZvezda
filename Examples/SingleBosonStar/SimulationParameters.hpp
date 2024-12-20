/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "BosonStarParams.hpp"
#include "ComplexPotential.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // Read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        pout() << "---------------------------------" << endl;

        // Gravitional constant
        pp.load("G_Newton", G_Newton, 1.0);

        // ######################################
        //  Single Boson Star Solver Parameters
        // ######################################

        // Boson Star initial data params
        pp.load("central_amplitude_CSF",
                bosonstar_params.central_amplitude_CSF);
        pp.load("phase", bosonstar_params.phase, 0.0);
        pp.load("gridpoints", bosonstar_params.gridpoints, 1000000);
        pp.load("BS_solver_psc", bosonstar_params.PSC, 2.0);
        pp.load("BS_solver_omc", bosonstar_params.OMC, 0.5);
        pp.load("BS_solver_verbosity", bosonstar_params.BS_solver_verbosity,
                false);
        pp.load("BS_solver_niter", bosonstar_params.niter, 17);

        pp.load("star_centre", bosonstar_params.star_centre, center);

        // Potential params
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
        pp.load("phi4_coeff", potential_params.phi4_coeff, 0.0);
        pp.load("solitonic", potential_params.solitonic, false);
        pp.load("sigma_solitonic", potential_params.sigma_solitonic, 0.2);

        positionA[0] = bosonstar_params.star_centre[0];
        positionA[1] = bosonstar_params.star_centre[1];
        positionA[2] = bosonstar_params.star_centre[2];

        pout() << "Star A is at x-position " << positionA[0] << endl;
        pout() << "Star A is at y-position " << positionA[1] << endl;
        pout() << "Star A is at z-position " << positionA[2] << endl;

#ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess,
                0.5 * bosonstar_params.mass);
#endif

        // Tagging
        pp.load("regrid_threshold_phi", regrid_threshold_phi);
        pp.load("regrid_threshold_chi", regrid_threshold_chi);

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("mass_write_extraction",
                mass_extraction_params.write_extraction, false);
        pp.load("num_mass_extraction_radii",
                mass_extraction_params.num_extraction_radii, 1);
        pp.load("mass_extraction_levels",
                mass_extraction_params.extraction_levels,
                mass_extraction_params.num_extraction_radii, 0);
        pp.load("mass_extraction_radii",
                mass_extraction_params.extraction_radii,
                mass_extraction_params.num_extraction_radii, 0.1);
        pp.load("num_points_phi_mass", mass_extraction_params.num_points_phi,
                2);
        pp.load("num_points_theta_mass",
                mass_extraction_params.num_points_theta, 4);
        pp.load("mass_extraction_center",
                mass_extraction_params.extraction_center,
                {0.5 * L, 0.5 * L, 0.5 * L});

        // Weyl extraction
        pp.load("activate_gw_extraction", activate_weyl_extraction, 0);

        // Work out the minimum extraction level
        auto min_extraction_level_it =
            mass_extraction_params.min_extraction_level();

        // Do we want to calculate and write the Noether Charge to a file
        pp.load("calculate_noether_charge", calculate_noether_charge, false);
    }

    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi;

    // Initial data for matter and potential
    double G_Newton;

    BosonStar_params_t bosonstar_params;
    ComplexPotential::params_t potential_params;

    // Mass extraction
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;

    int activate_weyl_extraction;

    // Do we want to write the Noether Charge to a file
    bool calculate_noether_charge;

    std::array<double, CH_SPACEDIM> positionA;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
