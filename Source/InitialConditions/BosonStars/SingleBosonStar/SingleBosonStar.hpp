/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SINGLEBOSONSTAR_HPP_
#define SINGLEBOSONSTAR_HPP_

#include "BosonStarHelperFunction.hpp"
#include "BosonStarParams.hpp"
#include "BosonStarSolver.hpp"
#include "Cell.hpp"
#include "ComplexPotential.hpp"
#include "ComplexScalarField.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "parstream.H"
#include "simd.hpp"

/*
 * This class constructs binary initial data for a single BS
 */

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

    BosonStarSolver m_1d_sol;

    double central_amplitude1, mass1, frequency1, radius1, compactness1;

  protected:
    double m_dx;
    BosonStar_params_t m_params_BosonStar;
    ComplexPotential::params_t m_params_potential;
};

#include "SingleBosonStar.impl.hpp"

#endif /* SINGLEBOSONSTAR_HPP_ */
