/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SINGLEBOSONSTAR_HPP_)
#error "This file should only be included through SingleBosonStar.hpp"
#endif

#ifndef SINGLEBOSONSTAR_IMPL_HPP_
#define SINGLEBOSONSTAR_IMPL_HPP_

inline SingleBosonStar::SingleBosonStar(
    BosonStar_params_t a_params_BosonStar,
    ComplexPotential::params_t a_params_potential, double a_dx)
    : m_dx(a_dx), m_params_BosonStar(a_params_BosonStar),
      m_params_potential(a_params_potential)
{
}

void SingleBosonStar::compute_1d_solution(const double max_r)
{
    try
    {
        // Set initial parameters for one star
        m_1d_sol.set_initialcondition_params(m_params_BosonStar,
                                             m_params_potential, max_r);
        m_1d_sol.main();
        central_amplitude1 = m_1d_sol.A[0];
        mass1 = m_1d_sol.boson_mass[m_params_BosonStar.gridpoints - 1];
        frequency1 = m_1d_sol.get_BSfrequency();
        radius1 = m_1d_sol.radius;
        compactness1 = m_1d_sol.compactness_value;
    }
    catch (std::exception &exception)
    {
        MayDay::Error("Yikes! I cannot compute the 1d solution in SingleBosonStar.impl.hpp file.");
    }
}

template <class data_t>
void SingleBosonStar::compute(Cell<data_t> current_cell) const
{
    using namespace TensorAlgebra;
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    current_cell.load_vars(vars);

    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_BosonStar.star_centre);

    // Import BS parameters
    double rapidity = m_params_BosonStar.BS_rapidity;
    bool antiboson = m_params_BosonStar.antiboson;
    double phase_offset = m_params_BosonStar.phase;

    BosonStarHelperFunction helper;

    // Compute star 1 and star 2 variables
    BosonStarHelperFunction::BS_3d_vars star1_vars = helper.compute_star_vars(coords, rapidity, 0., 0., m_1d_sol, phase_offset, antiboson);

    // Lapse and shift
    vars.lapse = star1_vars.lapse;
    vars.shift[0] = star1_vars.shiftx;

    // Scalar vars
    vars.phi_Re = star1_vars.phi_Re;
    vars.phi_Im = star1_vars.phi_Im;
    vars.Pi_Re = star1_vars.Pi_Re;
    vars.Pi_Im = star1_vars.Pi_Im;

    // Compute inverse metrics for star 1 and star 2
    auto gammaUU = compute_inverse_sym(star1_vars.gLL);

    // Define initial conformal factor
    vars.chi = pow(star1_vars.gLL[0][0] * star1_vars.gLL[1][1] * star1_vars.gLL[2][2], -1. / 3.);

    // Define initial trace of K and A_ij
    double one_third = 1. / 3.;
    FOR2(i, j) vars.h[i][j] = vars.chi * star1_vars.gLL[i][j];
    FOR2(i, j) vars.K = star1_vars.KLL[i][j] * gammaUU[i][j];
    FOR2(i, j)
    vars.A[i][j] = vars.chi * (star1_vars.KLL[i][j] - one_third * vars.K * star1_vars.gLL[i][j]);

    current_cell.store_vars(vars);
}

#endif /* SINGLEBOSONSTAR_IMPL_HPP_ */
