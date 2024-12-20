/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "SingleBosonStarLevel.hpp"
#include "BoxLoops.hpp"
#include "GammaCalculator.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "NewConstraints.hpp"
#include "NewMatterConstraints.hpp"

// For tag cells
#include "ComplexPhiAndChiExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "ComplexPotential.hpp"
#include "ComplexScalarField.hpp"
#include "ComputePack.hpp"
#include "SetValue.hpp"
#include "SingleBosonStar.hpp"

// For mass extraction
#include "ADMMass.hpp"
#include "ADMMassExtraction.hpp"
#include "EMTensor.hpp"
#include "MomFluxCalc.hpp"

// For GW extraction
#include "MatterWeyl4.hpp"
#include "WeylExtraction.hpp"

// For Noether Charge calculation
#include "NoetherCharge.hpp"
#include "SmallDataIO.hpp"

// for Chombo grid functions
#include "AMRReductions.hpp"

// Things to do at each advance step, after the RK4 is calculated
void SingleBosonStarLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance"),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void SingleBosonStarLevel::initialData()
{
    CH_TIME("SingleBosonStarLevel::initialData");
    if (m_verbosity)
        pout() << "SingleBosonStarLevel::initialData " << m_level << endl;

    // First initalise a BosonStar object
    SingleBosonStar boson_star(m_p.bosonstar_params, m_p.potential_params,
                               m_dx);

    // Initiate solver for 1D BS solutions
    boson_star.compute_1d_solution(4. * m_p.L);

    if (m_level == 0)
    {
        pout() << "Star 1 has A[0] " << boson_star.central_amplitude1
               << " mass " << boson_star.mass1 << " radius "
               << boson_star.radius1 << " and compactness "
               << boson_star.compactness1 << endl;
    }

    // First set everything to zero, as we do not want undefined values in
    // constraints. Then set initial conditions for a BS.
    BoxLoops::loop(make_compute_pack(SetValue(0.0), boson_star), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS, disable_simd());
}

// Things to do before outputting a checkpoint file
void SingleBosonStarLevel::preCheckpointLevel()
{
    CH_TIME("SingleBosonStarLevel::preCheckpointLevel");

    fillAllGhosts();
    ComplexPotential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(
        make_compute_pack(MatterConstraints<ComplexScalarFieldWithPotential>(
                              complex_scalar_field, m_dx, m_p.G_Newton, c_Ham,
                              Interval(c_Mom1, c_Mom3)),
                          NoetherCharge()),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do before outputting a plot file
void SingleBosonStarLevel::prePlotLevel()
{
    CH_TIME("SingleBosonStarLevel::prePlotLevel");

    fillAllGhosts();
    ComplexPotential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(
        make_compute_pack(MatterConstraints<ComplexScalarFieldWithPotential>(
                              complex_scalar_field, m_dx, m_p.G_Newton, c_Ham,
                              Interval(c_Mom1, c_Mom3)),
                          NoetherCharge()),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void SingleBosonStarLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                     const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ComplexScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these
    ComplexPotential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    MatterCCZ4RHS<ComplexScalarFieldWithPotential> my_ccz4_matter(
        complex_scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
        m_p.G_Newton);
    SetValue set_analysis_vars_zero(0.0, Interval(c_Pi_Im + 1, NUM_VARS - 1));
    auto compute_pack =
        make_compute_pack(my_ccz4_matter, set_analysis_vars_zero);
    BoxLoops::loop(compute_pack, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// Things to do at ODE update, after soln + rhs
void SingleBosonStarLevel::specificUpdateODE(GRLevelData &a_soln,
                                       const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void SingleBosonStarLevel::specificPostTimeStep()
{
    CH_TIME("SingleBosonStarLevel::specificPostTimeStep");

    bool first_step = (m_time == 0.0);

    // First compute diagnostics values on the grid
    fillAllGhosts();
    ComplexPotential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(MatterConstraints<ComplexScalarFieldWithPotential>(
                       complex_scalar_field, m_dx, m_p.G_Newton, c_Ham,
                       Interval(c_Mom1, c_Mom3)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    if (m_p.activate_mass_extraction == 1 &&
        m_level == m_p.mass_extraction_params.min_extraction_level())
    {
        if (m_verbosity)
        {
            pout() << "SingleBosonStarLevel::specificPostTimeStep:"
                      " Extracting mass."
                   << endl;
        }

        // Now refresh the interpolator and do the interpolation
        m_gr_amr.m_interpolator->refresh();
        ADMMassExtraction mass_extraction(m_p.mass_extraction_params, m_dt,
                                          m_time, first_step, m_restart_time);
        mass_extraction.execute_query(m_gr_amr.m_interpolator);
    }

    // Noether charge, max mod phi, min chi, constraint violations
    if (at_level_timestep_multiple(0))
    {
        BoxLoops::loop(NoetherCharge(), m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);
    }
    if (m_level == 0)
    {
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
        AMRReductions<VariableType::evolution> amr_reductions_ev(m_gr_amr);
        if (m_p.calculate_noether_charge)
        {
            // Compute volume weighted noether charge integral

            double noether_charge = amr_reductions.sum(c_N);
            SmallDataIO noether_charge_file("NoetherCharge", m_dt, m_time,
                                            m_restart_time, SmallDataIO::APPEND,
                                            first_step);
            noether_charge_file.remove_duplicate_time_data();
            if (m_time == 0.)
            {
                noether_charge_file.write_header_line({"Noether Charge"});
            }
            noether_charge_file.write_time_data_line({noether_charge});
        }

        // Compute the maximum of mod_phi and write it to a file
        double mod_phi_max = amr_reductions.max(c_mod_phi);
        SmallDataIO mod_phi_max_file("mod_phi_max", m_dt, m_time,
                                     m_restart_time, SmallDataIO::APPEND,
                                     first_step);
        mod_phi_max_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            mod_phi_max_file.write_header_line({"max mod phi"});
        }
        mod_phi_max_file.write_time_data_line({mod_phi_max});

        // Compute the min of \chi
        double min_chi = amr_reductions_ev.min(c_chi);
        SmallDataIO min_chi_file("min_chi", m_dt, m_time, m_restart_time,
                                 SmallDataIO::APPEND, first_step);
        min_chi_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            min_chi_file.write_header_line({"min chi"});
        }
        min_chi_file.write_time_data_line({min_chi});

        // Compute constraints 
        double L2_Ham = amr_reductions.norm(c_Ham, 2, true);
        double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3), 2, true);
        double L1_Ham = amr_reductions.norm(c_Ham, 1, true);
        double L1_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3), 1, true);
        SmallDataIO constraints_file("constraint_norms", m_dt, m_time,
                                     m_restart_time, SmallDataIO::APPEND,
                                     first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
        {
            constraints_file.write_header_line({
                "L^2_Ham",
                "L^2_Mom",
                "L^1_Ham",
                "L^1_Mom",
            });
        }
        constraints_file.write_time_data_line({L2_Ham, L2_Mom, L1_Ham, L1_Mom});
    }

#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        m_st_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
    }
#endif
}

void SingleBosonStarLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(ComplexPhiAndChiExtractionTaggingCriterion(m_dx, m_level,
                   m_p.extraction_params, m_p.regrid_threshold_phi,
                   m_p.regrid_threshold_chi, m_p.activate_extraction), current_state, tagging_criterion);
}
