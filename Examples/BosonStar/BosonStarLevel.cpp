/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "BosonStarLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"
#include "GammaCalculator.hpp"

// For RHS update
#include "MatterCCZ4.hpp"
#include "IntegratedMovingPunctureGauge.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"
#include "NewConstraints.hpp"

// For tag cells
#include "ComplexPhiAndChiExtractionTaggingCriterion.hpp"
#include "ChiandRhoTaggingCriterion.hpp"
#include "BosonChiPunctureExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "ComplexPotential.hpp"
#include "BosonStar.hpp"
#include "ComplexScalarField.hpp"
#include "SetValue.hpp"

// For mass extraction
#include "ADMMass.hpp"
//#include "Density.hpp"
#include "EMTensor.hpp"
#include "MomFluxCalc.hpp"
#include "SourceIntPreconditioner.hpp"
#include "ADMMassExtraction.hpp"

// For GW extraction
#include "MatterWeyl4.hpp"
#include "WeylExtraction.hpp"

// For Noether Charge calculation
#include "SmallDataIO.hpp"
#include "NoetherCharge.hpp"

// For Ang Mom Integrating
#include "AngMomFlux.hpp"

#include "ComputeWeightFunction.hpp"

// for chombo grid Functions
#include "AMRReductions.hpp"

// Things to do at each advance step, after the RK4 is calculated
void BosonStarLevel::specificAdvance()
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
void BosonStarLevel::initialData()
{
    CH_TIME("BosonStarLevel::initialData");
    if (m_verbosity)
        pout() << "BosonStarLevel::initialData " << m_level << endl;

    // First initalise a BosonStar object
    BosonStar boson_star(m_p.bosonstar_params, m_p.bosonstar2_params, m_p.potential_params,
                         m_p.G_Newton, m_dx, m_p.identical, m_verbosity);


    // the max radius the code might need to calculate out to is L*sqrt(3)
    boson_star.compute_1d_solution(4.*m_p.L);

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then  initial conditions for Boson Star
    BoxLoops::loop(make_compute_pack(SetValue(0.0), boson_star),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());
 
    BoxLoops::loop(GammaCalculator(m_dx),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS,
                   disable_simd());

    BoxLoops::loop(ComputeWeightFunction(m_p.bosonstar_params, m_p.bosonstar2_params, m_dx), m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS, disable_simd());

    fillAllGhosts();
    BoxLoops::loop(IntegratedMovingPunctureGauge(m_p.ccz4_params),
                 m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
    
}

// Things to do before outputting a checkpoint file
void BosonStarLevel::preCheckpointLevel()
{
    CH_TIME("BosonStarLevel::preCheckpointLevel");

     fillAllGhosts();
     Potential potential(m_p.potential_params);
     ComplexScalarFieldWithPotential complex_scalar_field(potential);
     BoxLoops::loop(make_compute_pack(
                     MatterWeyl4<ComplexScalarFieldWithPotential>(
                     complex_scalar_field,m_p.extraction_params.extraction_center,
                     m_dx, m_p.formulation, m_p.G_Newton),
                     MatterConstraints<ComplexScalarFieldWithPotential>(
                     complex_scalar_field, m_dx, m_p.G_Newton, c_Ham,
                     Interval(c_Mom1, c_Mom3)), NoetherCharge(),
                     EMTensor<ComplexScalarFieldWithPotential>(
                     complex_scalar_field, m_dx, c_rho, Interval(c_s1,c_s3),
                     Interval(c_s11,c_s33))),
                     m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

}

// Things to do before outputting a plot file
void BosonStarLevel::prePlotLevel()
{
    CH_TIME("BosonStarLevel::prePlotLevel");

      fillAllGhosts();
      Potential potential(m_p.potential_params);
      ComplexScalarFieldWithPotential complex_scalar_field(potential);
      BoxLoops::loop(make_compute_pack(
                      MatterWeyl4<ComplexScalarFieldWithPotential>(
                      complex_scalar_field,m_p.extraction_params.extraction_center,
                      m_dx, m_p.formulation, m_p.G_Newton),
                      MatterConstraints<ComplexScalarFieldWithPotential>(
                      complex_scalar_field, m_dx, m_p.G_Newton, c_Ham,
                      Interval(c_Mom1, c_Mom3)), NoetherCharge(),
                      EMTensor<ComplexScalarFieldWithPotential>(
                      complex_scalar_field, m_dx, c_rho, Interval(c_s1,c_s3),
                      Interval(c_s11,c_s33))),
                      m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

}

// Things to do in RHS update, at each RK4 step
void BosonStarLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()), a_soln,
        a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ComplexScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    //MatterCCZ4RHS<ComplexScalarFieldWithPotential> my_ccz4_matter(
    //    complex_scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
    //    m_p.G_Newton);
    MatterCCZ4RHS<ComplexScalarFieldWithPotential, IntegratedMovingPunctureGauge, FourthOrderDerivatives> my_ccz4_matter(
         complex_scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
         m_p.G_Newton);
    SetValue set_analysis_vars_zero(0.0, Interval(c_Pi_Im + 1, NUM_VARS - 1));
    auto compute_pack =
        make_compute_pack(my_ccz4_matter, set_analysis_vars_zero);
    BoxLoops::loop(compute_pack, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// Things to do at ODE update, after soln + rhs
void BosonStarLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void BosonStarLevel::specificPostTimeStep()
{
    CH_TIME("BosonStarLevel::specificPostTimeStep");

    bool first_step = (m_time == 0.0);

    // First compute the ADM Mass integrand values on the grid
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    auto weyl4_adm_compute_pack =
               make_compute_pack(MatterWeyl4<ComplexScalarFieldWithPotential>(
               complex_scalar_field,m_p.extraction_params.extraction_center,
               m_dx, m_p.formulation, m_p.G_Newton), ADMMass(m_p.center, m_dx));
    BoxLoops::loop(weyl4_adm_compute_pack, m_state_new, m_state_diagnostics,
                        EXCLUDE_GHOST_CELLS);
    BoxLoops::loop(MatterConstraints<ComplexScalarFieldWithPotential>(
                        complex_scalar_field, m_dx, m_p.G_Newton, c_Ham,
                        Interval(c_Mom1, c_Mom3)), m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    
    if (m_p.activate_weyl_extraction == 1 &&
       at_level_timestep_multiple(m_p.extraction_params.min_extraction_level()))
    {
        CH_TIME("BosonStarLevel::doAnalysis::Weyl4&ADMMass");
        
        // Do the extraction on the min extraction level
        if (m_level == m_p.extraction_params.min_extraction_level())
        {
            if (m_verbosity)
            {
                pout() << "BinaryBSLevel::specificPostTimeStep:"
                          " Extracting gravitational waves." << endl;
            }


            // Refresh the interpolator and do the interpolation
            m_gr_amr.m_interpolator->refresh();
            WeylExtraction gw_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gw_extraction.execute_query(m_gr_amr.m_interpolator);
        }
    }

    if (m_p.activate_mass_extraction == 1 &&
        m_level == m_p.mass_extraction_params.min_extraction_level())
    {
        if (m_verbosity)
        {
            pout() << "BinaryBSLevel::specificPostTimeStep:"
                      " Extracting mass." << endl;
        }

        // Now refresh the interpolator and do the interpolation
        m_gr_amr.m_interpolator->refresh();
        ADMMassExtraction mass_extraction(m_p.mass_extraction_params, m_dt,
                                    m_time, first_step,
                                    m_restart_time);
        mass_extraction.execute_query(m_gr_amr.m_interpolator);
    }

    // noether charge, max mod phi, min chi, constraint violations
    if (at_level_timestep_multiple(0))
    {
        BoxLoops::loop(NoetherCharge(), m_state_new, m_state_diagnostics,
                  EXCLUDE_GHOST_CELLS);
    }
    if (m_level == 0)
    {
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
        if (m_p.calculate_noether_charge)
        {
            //noether charge should be calculated pre-check and pre plot
            //so automatically here

            // compute integrated volume weighted noether charge integral

            double noether_charge = amr_reductions.sum(c_N);
            SmallDataIO noether_charge_file("NoetherCharge", m_dt, m_time,
                                            m_restart_time,
                                            SmallDataIO::APPEND,
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
                                     m_restart_time,
                                     SmallDataIO::APPEND,
                                     first_step);
        mod_phi_max_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            mod_phi_max_file.write_header_line({"max mod phi"});
        }
        mod_phi_max_file.write_time_data_line({mod_phi_max});


        // Compute the min of chi and write it to a file
        double min_chi = amr_reductions.min(c_chi);
        SmallDataIO min_chi_file("min_chi", m_dt, m_time,
                                     m_restart_time,
                                     SmallDataIO::APPEND,
                                     first_step);
        min_chi_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            min_chi_file.write_header_line({"min chi"});
        }
        min_chi_file.write_time_data_line({min_chi});


        // constraints calculated pre check and pre plot so done here already

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
            constraints_file.write_header_line({"L^2_Ham", "L^2_Mom", "L^1_Ham", "L^1_Mom",});
        }
        constraints_file.write_time_data_line({L2_Ham, L2_Mom, L1_Ham, L1_Mom});
    }

    if (m_p.do_star_track && m_level == m_p.star_track_level)
    {
	pout() << "Running a star tracker now" << endl;
        // if at restart time read data from dat file,
        // will default to param file if restart time is 0
        if (fabs(m_time - m_restart_time) < m_dt * 1.1)
        {
            m_st_amr.m_star_tracker.read_old_centre_from_dat(
                "StarCentres", m_dt, m_time, m_restart_time, first_step);
        }
        m_st_amr.m_star_tracker.update_star_centres(m_dt);
        m_st_amr.m_star_tracker.write_to_dat("StarCentres", m_dt, m_time,
                                             m_restart_time, first_step);
    }

    //if (m_p.do_flux_integration && m_level==m_p.angmomflux_params.extraction_level)
    double temp_dx;
    // if (m_p.do_flux_integration && at_level_timestep_multiple(m_p.flux_extraction_level))
    // {
    //     CH_TIME("BosonStarLevel::doAnalysis::FphiSphi");
    //     Potential potential(m_p.potential_params);
    //     ComplexScalarFieldWithPotential complex_scalar_field(potential);

    //     // execute only on the finest level
    //     if (m_level==m_p.flux_extraction_level)
    //     //if (m_level==m_p.max_level)
    //     {
    //         double S_phi_integral; // integral of angmomsource
    //         double Q_phi_integral; // integral of angmomsource
    //         std::vector<AMRLevel *> all_level_ptrs = getAMRLevelHierarchy().stdVector();
    //         std::vector<double> S_phi_integrals(m_p.angmomflux_params.num_extraction_radii); // vector storing all integrals
    //         std::vector<double> Q_phi_integrals(m_p.angmomflux_params.num_extraction_radii); // vector storing all integrals

    //         //temp_dx = m_p.coarsest_dx;
    //         // fill grid with angmom variables on every level
    //         for (auto level_ptr : all_level_ptrs)
    //         {
    //             BosonStarLevel *bs_level_ptr =
    //                 dynamic_cast<BosonStarLevel *>(level_ptr);
    //             if (bs_level_ptr == nullptr)
    //             {
    //                 break;
    //             }
    //             temp_dx = bs_level_ptr->m_dx;
    //             BoxLoops::loop(EMTensor_and_mom_flux<ComplexScalarFieldWithPotential>(
    //             complex_scalar_field, temp_dx, m_p.L, m_p.angmomflux_params.center,
    //                                  c_Fphi_flux, c_Sphi_source, c_Qphi_density,
    //                                                  c_rho, Interval(c_s1,c_s3),
    //                          Interval(c_s11,c_s33)),  bs_level_ptr->m_state_new,
    //                             bs_level_ptr->m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    //         }
    //         // loop through levels and fill ghosts
    //         for (auto level_ptr : all_level_ptrs)
    //         {
    //             BosonStarLevel *bs_level_ptr =
    //                 dynamic_cast<BosonStarLevel *>(level_ptr);
    //             if (bs_level_ptr == nullptr)
    //             {
    //                 break;
    //             }
    //             bs_level_ptr->fillAllGhosts();
    //         }

            // for (int i=m_p.angmomflux_params.num_extraction_radii-1; i>=0; i--)
            // {
            //     temp_dx = m_p.coarsest_dx;
            //     for (auto level_ptr : all_level_ptrs)
            //     {
            //         BosonStarLevel *bs_level_ptr =
            //             dynamic_cast<BosonStarLevel *>(level_ptr);
            //         if (bs_level_ptr == nullptr)
            //         {
            //             break;
            //         }

            //       // set angmomsource and density to zero outside of extraction radii
            //       temp_dx = bs_level_ptr->m_dx;
            //       BoxLoops::loop(SourceIntPreconditioner<ComplexScalarFieldWithPotential>
            //   (complex_scalar_field, temp_dx, m_p.L, m_p.angmomflux_params.center,
            //                                      c_Sphi_source, c_Qphi_density,
            //                         m_p.angmomflux_params.extraction_radii[i]),
            //               bs_level_ptr->m_state_new, bs_level_ptr->m_state_diagnostics,
            //                                               INCLUDE_GHOST_CELLS);
            //     }
            //     // old code passed this coarse dx, maybe this doesnt work now?
            //     AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            //     S_phi_integral = amr_reductions.sum(c_Sphi_source);
            //     S_phi_integrals[i] = S_phi_integral;
            //     Q_phi_integral = amr_reductions.sum(c_Qphi_density);
            //     Q_phi_integrals[i] = Q_phi_integral;
            // }


            // save the Source integral to dat file
            // std::vector<string> title_line(m_p.angmomflux_params.num_extraction_radii);
            // string dummy_string;
            // for (int j=0; j<m_p.angmomflux_params.num_extraction_radii; j++)
            // {
            //     dummy_string = "r = " + to_string(m_p.angmomflux_params.extraction_radii[j]);
            //     title_line[j] = dummy_string;
            // }

            // SmallDataIO angmomsource_file("AngMomSource", m_dt, m_time,
            //                               m_restart_time,
            //                               SmallDataIO::APPEND,
            //                               first_step);

            // if (m_time > 0) angmomsource_file.remove_duplicate_time_data();

            // if (m_time == 0.)
            // {
            //     angmomsource_file.write_header_line(title_line);
            // }

            // angmomsource_file.write_time_data_line(S_phi_integrals);


            // // save the Density integral to dat file
            // std::vector<string> title_line2(m_p.angmomflux_params.num_extraction_radii);
            // string dummy_string2;
            // for (int j=0; j<m_p.angmomflux_params.num_extraction_radii; j++)
            // {
            //     dummy_string2 = "r = " + to_string(m_p.angmomflux_params.extraction_radii[j]);
            //     title_line2[j] = dummy_string2;
            // }

    //         SmallDataIO density_file("AngMomDensity", m_dt, m_time,
    //                                       m_restart_time,
    //                                       SmallDataIO::APPEND,
    //                                       first_step);

    //         if (m_time > 0) density_file.remove_duplicate_time_data();

    //         if (m_time == 0.)
    //         {
    //             density_file.write_header_line(title_line2);
    //         }

    //         density_file.write_time_data_line(Q_phi_integrals);



    //         // Refresh the interpolator and do the interpolation
    //         m_gr_amr.m_interpolator->refresh();
    //         // setup and do angmomflux integral
    //         AngMomFlux ang_mom_flux(m_p.angmomflux_params,m_time,m_dt,m_restart_time,first_step);
    //         ang_mom_flux.run(m_gr_amr.m_interpolator);
    //     }
    // }
}

void BosonStarLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{

//BoxLoops::loop(ComplexPhiAndChiExtractionTaggingCriterion(m_dx, m_level,
//                   m_p.mass_extraction_params, m_p.regrid_threshold_phi,
//                   m_p.regrid_threshold_chi), current_state, tagging_criterion);

//    BoxLoops::loop(ChiandRhoTaggingCriterion(m_dx, m_level,
//                    m_p.mass_extraction_params, m_p.regrid_threshold_rho,
//                    m_p.regrid_threshold_chi), current_state, tagging_criterion);

// Be aware of the tagging here, you may want to change it, depending on your problem of interest. Below tagging for when the tracking is activated is 'intense' and specific to binary inspirals. 	
    
     if (m_p.do_star_track == true)
    {
        const vector<double> puncture_radii = {m_p.tag_radius_A,
                                                m_p.tag_radius_B};
        const vector<double> puncture_masses = {m_p.bosonstar_params.mass,
                                                m_p.bosonstar2_params.mass};

        const std::vector<double> star_coords =
            m_st_amr.m_star_tracker.get_puncture_coords();
	
        BoxLoops::loop(BosonChiPunctureExtractionTaggingCriterion(
                           m_dx, m_level, m_p.tag_horizons_max_levels,
                       m_p.tag_punctures_max_levels, m_p.extraction_params,
                           star_coords, m_p.activate_extraction,
                           m_p.do_star_track, puncture_radii, puncture_masses, m_p.tag_buffer),
                      current_state, tagging_criterion);	
   }
    else
   {
        BoxLoops::loop(ChiandRhoTaggingCriterion(m_dx, m_level,
                   m_p.mass_extraction_params, m_p.regrid_threshold_rho,
                   m_p.regrid_threshold_chi), current_state, tagging_criterion);
   }

}
