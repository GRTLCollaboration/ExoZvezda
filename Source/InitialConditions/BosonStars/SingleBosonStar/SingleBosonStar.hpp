/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SINGLEBOSONSTAR_HPP_
#define SINGLEBOSONSTAR_HPP_

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
#include "parstream.H" //gives pout
#include "simd.hpp"

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class SingleBosonStar
{

  public:
    //! The constructor
    SingleBosonStar(BosonStar_params_t a_params_BosonStar,
                    ComplexPotential::params_t a_params_potential, double a_dx);

    //! Computes the 1d solution and stores in m_1d_sol
    void compute_1d_solution(const double max_r);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    BosonStarSolution m_1d_sol;

    double central_amplitude1, mass1, frequency1, radius1, compactness1;
    
  protected:
    double m_dx;
    BosonStar_params_t m_params_BosonStar;  //!< The complex scalar field params
    ComplexPotential::params_t m_params_potential; //!< The potential params
};

#include "SingleBosonStar.impl.hpp"

#endif /* SINGLEBOSONSTAR_HPP_ */
