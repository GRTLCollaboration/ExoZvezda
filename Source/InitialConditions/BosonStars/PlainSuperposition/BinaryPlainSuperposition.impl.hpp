/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BINARYPLAINSUPERPOSITION_HPP_)
#error
"This file should only be included through BinaryPlainSuperposition.impl.hpp"
#endif

#ifndef BINARYPLAINSUPERPOSITION_IMPL_HPP_
#define BINARYPLAINSUPERPOSITION_IMPL_HPP_

    inline BinaryPlainSuperposition::BinaryPlainSuperposition(
        BosonStar_params_t a_params_BosonStar,
        BosonStar_params_t a_params_BosonStar2,
        ComplexPotential::params_t a_params_potential, double a_dx)
    : m_dx(a_dx), m_params_BosonStar(a_params_BosonStar),
      m_params_BosonStar2(a_params_BosonStar2),
      m_params_potential(a_params_potential)
{
}

// Solution is integrated up to large radius given by max_r.
void BinaryPlainSuperposition::compute_1d_solution(const double max_r)
{
    try
    {
        // Compute for the 1st star
        m_1d_sol.set_initialcondition_params(m_params_BosonStar,
                                             m_params_potential, max_r);
        m_1d_sol.main();
        central_amplitude1 = m_1d_sol.A[0];
        mass1 = m_1d_sol.boson_mass[m_params_BosonStar.gridpoints - 1];
        frequency1 = m_1d_sol.get_BSfrequency();
        radius1 = m_1d_sol.radius;
        compactness1 = m_1d_sol.compactness_value;

        // Compute for the 2nd star
        m_1d_sol2.set_initialcondition_params(m_params_BosonStar2,
                                              m_params_potential, max_r);
        m_1d_sol2.main();
        central_amplitude2 = m_1d_sol2.A[0];
        mass2 = m_1d_sol2.boson_mass[m_params_BosonStar.gridpoints - 1];
        frequency2 = m_1d_sol2.get_BSfrequency();
        radius2 = m_1d_sol2.radius;
        compactness2 = m_1d_sol2.compactness_value;
    }
    catch (std::exception &exception)
    {
        MayDay::Error("Yikes! I cannot compute the 1d solution in "
                      "BinaryPlainSuperposition.impl.hpp file.");
    }
}

template <class data_t>
void BinaryPlainSuperposition::compute(Cell<data_t> current_cell) const
{
    using namespace TensorAlgebra;
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    current_cell.load_vars(vars);

    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_BosonStar.star_centre);

    // Import BS parameters and option of whether this is a BS binary or BS-BH
    // binary
    double rapidity = m_params_BosonStar.BS_rapidity;
    double rapidity2 = m_params_BosonStar2.BS_rapidity;
    bool antiboson = m_params_BosonStar.antiboson;
    double phase_offset = m_params_BosonStar.phase;
    double separation = m_params_BosonStar.BS_separation;
    double impact_parameter = m_params_BosonStar.BS_impact_parameter;
    double q = m_params_BosonStar.mass_ratio;

    BosonStarHelperFunction helper;

    // Compute star 1 and star 2 variables
    BosonStarHelperFunction::BS_3d_vars star1_vars = helper.compute_star_vars(
        coords, rapidity, -q * separation / (q + 1.),
        q * impact_parameter / (q + 1.), m_1d_sol, phase_offset, false);
    BosonStarHelperFunction::BS_3d_vars star2_vars = helper.compute_star_vars(
        coords, -rapidity2, separation / (q + 1.), -impact_parameter / (q + 1.),
        m_1d_sol2, 0.0, antiboson);

    // Superpose shift
    vars.shift[0] += star1_vars.shiftx + star2_vars.shiftx;

    // Superpose scalar vars
    vars.phi_Re += star1_vars.phi_Re + star2_vars.phi_Re;
    vars.phi_Im += star1_vars.phi_Im + star2_vars.phi_Im;
    vars.Pi_Re += star1_vars.Pi_Re + star2_vars.Pi_Re;
    vars.Pi_Im += star1_vars.Pi_Im + star2_vars.Pi_Im;

    // Compute inverse metrics for star 1 and star 2
    auto gammaUU_1 = compute_inverse_sym(star1_vars.gLL);
    auto gammaUU_2 = compute_inverse_sym(star2_vars.gLL);

    Tensor<2, data_t> gammaLL,
        KLL; // to store superposed metric and extrinsic curvature
    FOR(i, j)
    {
        gammaLL[i][j] = 0.0;
        KLL[i][j] = 0.0;
    }

    // Superpose metrics
    gammaLL[0][0] = star1_vars.gLL[0][0] + star2_vars.gLL[0][0] - 1.0;
    gammaLL[1][1] = star1_vars.gLL[1][1] + star2_vars.gLL[1][1] - 1.0;
    gammaLL[2][2] = star1_vars.gLL[2][2] + star2_vars.gLL[2][2] - 1.0;

    // Compute inverse
    auto gammaUU = compute_inverse_sym(gammaLL);

    // Define initial conformal factor
    vars.chi = pow(gammaLL[0][0] * gammaLL[1][1] * gammaLL[2][2], -1. / 3.);

    // Define initial lapse
    vars.lapse += sqrt(star1_vars.lapse * star1_vars.lapse +
                       star2_vars.lapse * star2_vars.lapse - 1.);

    // Define initial trace of K and A_ij
    double one_third = 1. / 3.;
    FOR2(i, j) vars.h[i][j] = vars.chi * gammaLL[i][j];
    FOR4(i, j, k, l)
    KLL[i][j] += gammaLL[i][l] * (gammaUU_1[l][k] * star1_vars.KLL[k][j] +
                                  gammaUU_2[l][k] * star2_vars.KLL[k][j]);
    FOR2(i, j) vars.K += KLL[i][j] * gammaUU[i][j];
    FOR2(i, j)
    vars.A[i][j] = vars.chi * (KLL[i][j] - one_third * vars.K * gammaLL[i][j]);

    current_cell.store_vars(vars);
}

#endif /* BINARYPLAINSUPERPOSITION_IMPL_HPP_ */
