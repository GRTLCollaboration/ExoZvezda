/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "CH_Timer.H"
#include "parstream.H" //Gives us pout()

// System includes
#include <chrono>
#include <iostream>

#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "BBSEqualMassFixLevel.hpp"

// Star tracking
#include "STAMR.hpp"

int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    // The line below selects the problem that is simulated
    // (To simulate a different problem, define a new child of AMRLevel
    // and an associated LevelFactory)
    STAMR st_amr;

    st_amr.m_star_tracker.initialise_star_tracking(
        sim_params.do_star_track, sim_params.number_of_stars,
        {sim_params.positionA, sim_params.positionB}, sim_params.star_points,
        sim_params.star_track_width_A, sim_params.star_track_width_B,
        sim_params.star_track_direction_of_motion);
    DefaultLevelFactory<BBSEqualMassFixLevel> boson_star_level_fact(st_amr,
                                                                    sim_params);
    setupAMRObject(st_amr, boson_star_level_fact);

    // Instantiate AMR interpolator for mass/GW extraction
    AMRInterpolator<Lagrange<4>> interpolator(
        st_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    st_amr.set_interpolator(&interpolator);

    // Add a scheduler to GRAMR which just calls doAnalysis on every AMRLevel
    // at time 0. It is called later in postTimeStep
    // RefCountedPtr<CallDoAnalysis> call_do_analysis_ptr(new CallDoAnalysis);
    // RefCountedPtr<Scheduler> scheduler_ptr(new Scheduler);
    // scheduler_ptr->schedule(call_do_analysis_ptr, sim_params.max_steps);
    // st_amr.schedule(scheduler_ptr);

#ifdef USE_AHFINDER
    if (sim_params.AH_activate)
    {
        AHSurfaceGeometry sph1(sim_params.positionA);
        AHSurfaceGeometry sph2(sim_params.positionB);

        st_amr.m_ah_finder.add_ah(sph1, sim_params.AH_1_initial_guess,
                                  sim_params.AH_params);
        st_amr.m_ah_finder.add_ah(sph2, sim_params.AH_2_initial_guess,
                                  sim_params.AH_params);
        st_amr.m_ah_finder.add_ah_merger(0, 1, sim_params.AH_params);
    }
#endif

    using Clock = std::chrono::steady_clock;
    using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;

    std::chrono::time_point<Clock> start_time = Clock::now();

    // Add a scheduler to call specificPostTimeStep on every AMRLevel at t=0
    auto task = [](GRAMRLevel *level)
    {
        if (level->time() == 0.)
            level->specificPostTimeStep();
    };
    // call 'now' really now
    MultiLevelTaskPtr<> call_task(task);
    call_task.execute(st_amr);

    // Engage! Run the evolution.
    st_amr.run(sim_params.stop_time, sim_params.max_steps);

    auto now = Clock::now();
    auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
    pout() << "Total simulation time (mins): " << duration.count() << ".\n";

    st_amr.conclude();

    // Write Chombo timer report
    CH_TIMER_REPORT();

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChombo(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
