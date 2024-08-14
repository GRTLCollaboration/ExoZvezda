/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARLEVEL_HPP_
#define BOSONSTARLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "ComplexPotential.hpp"
#include "ComplexScalarField.hpp"
#include "STAMR.hpp"

//!  A class for the evolution of a single boson star.
/*!
     The class takes some initial data for a boson star (i.e. a compact complex
     scalar field configuration) and evolves it using the CCZ4 equations.
*/
class BosonStarLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<BosonStarLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    STAMR &m_st_amr = dynamic_cast<STAMR &>(m_gr_amr);

    // Typedef for scalar field
    typedef ComplexScalarField<Potential> ComplexScalarFieldWithPotential;

    //! Things to do at the end of the advance step, after RK4 calculation
    virtual void specificAdvance() override;

    //! Initialize data for the field and metric variables
    virtual void initialData() override;

    //! routines to do before outputing checkpoint file
    virtual void preCheckpointLevel() override;

    //! routines to do before outputing plot file
    virtual void prePlotLevel() override;

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    //! Things to do in UpdateODE step, after soln + rhs update
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    //! Tell Chombo how to tag cells for regridding
    virtual void computeTaggingCriterion(
        FArrayBox &tagging_criterion, const FArrayBox &current_state,
        const FArrayBox &current_state_diagnostics) override;

    //! Things to do for analysis after each timestep and at the start
    virtual void specificPostTimeStep() override;
};

#endif /* BOSONSTARLEVEL_HPP_ */
