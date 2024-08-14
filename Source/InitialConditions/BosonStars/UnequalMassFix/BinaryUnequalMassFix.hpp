/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYUNEQUALMASSFIX_HPP_
#define BINARYUNEQUALMASSFIX_HPP_

#include "BosonStarParams.hpp"
#include "BosonStarSolution.hpp"
#include "Cell.hpp"
#include "ComplexPotential.hpp"
#include "ComplexScalarField.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "WeightFunction.hpp"
#include "parstream.H" //gives pout
#include "simd.hpp"

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class BinaryUnequalMassFix
{

  public:
    //! The constructor
    BinaryUnequalMassFix(BosonStar_params_t a_params_BosonStar,
                         BosonStar_params_t a_params_BosonStar2,
                         Potential::params_t a_params_potential, double a_dx);

    //! Computes the 1d solution and stores in m_1d_sol
    void compute_1d_solution(const double max_r);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    BosonStarSolution m_1d_sol;
    BosonStarSolution m_1d_sol2;

    // The object that stores the solution found by the 1d ODE integrator */

  protected:
    double m_dx;
    BosonStar_params_t m_params_BosonStar;
    BosonStar_params_t m_params_BosonStar2; //!< The complex scalar field params
    Potential::params_t m_params_potential; //!< The potential params
};

#include "BinaryUnequalMassFix.impl.hpp"

#endif /* BINARYUNEQUALMASSFIX_HPP_ */
