/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BINARYUNEQUALMASSFIX_HPP_)
#error "This file should only be included through BinaryUnequalMassFix.hpp"
#endif

#ifndef BINARYUNEQUALMASSFIX_IMPL_HPP_
#define BINARYUNEQUALMASSFIX_IMPL_HPP_

inline BinaryUnequalMassFix::BinaryUnequalMassFix(
    BosonStar_params_t a_params_BosonStar,
    BosonStar_params_t a_params_BosonStar2,
    ComplexPotential::params_t a_params_potential, double a_dx)
    : m_dx(a_dx), m_params_BosonStar(a_params_BosonStar),
      m_params_BosonStar2(a_params_BosonStar2),
      m_params_potential(a_params_potential)
{
}

void BinaryUnequalMassFix::compute_1d_solution(const double max_r)
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
                      "BinaryUnequalMassFix.impl.hpp file.");
    }
}

template <class data_t>
void BinaryUnequalMassFix::compute(Cell<data_t> current_cell) const
{
    using namespace TensorAlgebra;
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    current_cell.load_vars(vars);

    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_BosonStar.star_centre);

    // Import BS parameters
    double rapidity = m_params_BosonStar.BS_rapidity;
    double rapidity2 = m_params_BosonStar2.BS_rapidity;
    bool antiboson = m_params_BosonStar.antiboson;
    double phase_offset = m_params_BosonStar.phase;
    double separation = m_params_BosonStar.BS_separation;
    double impact_parameter = m_params_BosonStar.BS_impact_parameter;
    double q = m_params_BosonStar.mass_ratio;
    double radius_width1 = m_params_BosonStar.radius_width1;
    double radius_width2 = m_params_BosonStar.radius_width2;
    int conformal_power = m_params_BosonStar.conformal_factor_power;

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

    // Start with plain superposed metrics
    gammaLL[0][0] = star1_vars.gLL[0][0] + star2_vars.gLL[0][0] - 1.0;
    gammaLL[1][1] = star1_vars.gLL[1][1] + star2_vars.gLL[1][1] - 1.0;
    gammaLL[2][2] = star1_vars.gLL[2][2] + star2_vars.gLL[2][2] - 1.0;

    // Compute inverse
    auto gammaUU = compute_inverse_sym(gammaLL);

    // We are now ready to start applying the unequal-mass fix, but we need to
    // calculate a couple of terms
    Tensor<2, data_t> correctionLL;  // correction needed for star A
    Tensor<2, data_t> correctionLL2; // correction needed for star B
    Tensor<2, data_t> superposed_1;  // gamma_{ij}(x_A) from plain superposition
    Tensor<2, data_t> superposed_2;  // gamma_{ij}(x_B) from plain superposition

    FOR(i, j)
    {
        correctionLL[i][j] = 0.0;
        correctionLL2[i][j] = 0.0;
        superposed_1[i][j] = 0.0;
        superposed_2[i][j] = 0.0;
    }

    // Compute the effect of object 1 on object 2 and write into correctionLL
    double c_ = cosh(rapidity);
    double s_ = sinh(rapidity);
    double x_p = (-separation) * c_;
    double z_p = 0.;
    double y_p = impact_parameter;
    double r_p = sqrt(x_p * x_p + y_p * y_p + z_p * z_p);
    double omega_p = m_1d_sol.get_lapse_interp(r_p);
    double psi_p = m_1d_sol.get_psi_interp(r_p);
    double pc_os_p = psi_p * psi_p * c_ * c_ - omega_p * omega_p * s_ * s_;

    correctionLL[1][1] = psi_p * psi_p;
    correctionLL[2][2] = psi_p * psi_p;
    correctionLL[0][0] = pc_os_p;

    // Compute the effect of object 2 on object 1 and write into correctionLL2
    double c2_ = cosh(-rapidity2);
    double s2_ = sinh(-rapidity2);
    double x_p2 = (separation)*c2_;
    double z_p2 = 0.;
    double y_p2 = -impact_parameter;
    double r_p2 = sqrt(x_p2 * x_p2 + y_p2 * y_p2 + z_p2 * z_p2);
    double omega_p2 = m_1d_sol2.get_lapse_interp(r_p2);
    double psi_p2 = m_1d_sol2.get_psi_interp(r_p2);
    double pc_os_p2 =
        psi_p2 * psi_p2 * c2_ * c2_ - omega_p2 * omega_p2 * s2_ * s2_;

    correctionLL2[1][1] = psi_p2 * psi_p2;
    correctionLL2[2][2] = psi_p2 * psi_p2;
    correctionLL2[0][0] = pc_os_p2;

    // Compute the values of the metric of star 1 at its centre
    double r_11 = 0.;
    double omega_11 = m_1d_sol.get_lapse_interp(r_11);
    double psi_11 = m_1d_sol.get_psi_interp(r_11);
    double pc_os_11 = psi_11 * psi_11 * cosh(rapidity) * cosh(rapidity) -
                      omega_11 * omega_11 * sinh(rapidity) * sinh(rapidity);

    // metric components of \gamma_A(x_A)
    double g_zz_11 = psi_11 * psi_11;
    double g_yy_11 = psi_11 * psi_11;
    double g_xx_11 = pc_os_11;

    // Compute the values of the metric of star 2 at its centre
    double r_22 = 0.;
    double omega_22 = m_1d_sol2.get_lapse_interp(r_22);
    double psi_22 = m_1d_sol2.get_psi_interp(r_22);
    double pc_os_22 = psi_22 * psi_22 * cosh(-rapidity2) * cosh(-rapidity2) -
                      omega_22 * omega_22 * sinh(-rapidity2) * sinh(-rapidity2);

    // metric components of \gamma_B(x_B)
    double g_zz_22 = psi_22 * psi_22;
    double g_yy_22 = psi_22 * psi_22;
    double g_xx_22 = pc_os_22;

    // This is \gamma_{ij}(x_A) = \gamma_A(x_A) + \gamma_B(x_A) - 1
    superposed_1[0][0] = g_xx_11 + correctionLL2[0][0] - 1.;
    superposed_1[1][1] = g_yy_11 + correctionLL2[1][1] - 1.;
    superposed_1[2][2] = g_zz_11 + correctionLL2[2][2] - 1.;

    // This  is \gamma_{ij}(x_B) = \gamma_B(x_B) + \gamma_A(x_B) - 1
    superposed_2[0][0] = g_xx_22 + correctionLL[0][0] - 1.;
    superposed_2[1][1] = g_yy_22 + correctionLL[1][1] - 1.;
    superposed_2[2][2] = g_zz_22 + correctionLL[2][2] - 1.;

    // Eq. (A.4) of https://arxiv.org/abs/2212.08023
    double n_power = conformal_power / 12.0;
    // Conformal factor from plain superposition
    double chi_plain =
        pow(gammaLL[0][0] * gammaLL[1][1] * gammaLL[2][2], n_power);

    // Now, we correct the conformal factor using Eqs.(19)-(20) of
    // https://arxiv.org/abs/2212.08023 (in the text we use \lambda instead of
    // \chi)

    // This is \chi(x_A) of Eq.(19)
    double chi_1 = pow(
        superposed_1[0][0] * superposed_1[1][1] * superposed_1[2][2], n_power);
    // This is \chi(x_B) of Eq.(20)
    double chi_2 = pow(
        superposed_2[0][0] * superposed_2[1][1] * superposed_2[2][2], n_power);

    // This is \chi^A(x_A) of Eq.(19)
    double chi1_1 = pow(g_xx_11 * g_yy_11 * g_zz_11, n_power);
    // This is \chi^B(x_B) of Eq.(20)
    double chi2_2 = pow(g_xx_22 * g_yy_22 * g_zz_22, n_power);

    // This is \delta_A of Eq.(19)
    double delta_1 = chi1_1 - chi_1;
    // This is \delta_B of Eq.(20)
    double delta_2 = chi2_2 - chi_2;

    // Initialise weight function calculation
    WeightFunction weight;

    // Find all the profile functions needed
    double profile1 = weight.profile_chi(
        (coords.x - q * separation / (q + 1)) * cosh(rapidity),
        coords.y + q * impact_parameter / (q + 1.), coords.z,
        radius_width1); // w_A(x)
    double profile2 =
        weight.profile_chi((coords.x + separation / (q + 1)) * cosh(-rapidity2),
                           coords.y - impact_parameter / (q + 1.), coords.z,
                           radius_width2); // w_B(x)

    double profile_11 =
        weight.profile_chi(0., 0., 0., radius_width1); // w_A(x_A)
    double argument_xB_xA =
        (separation / (q + 1)) *
        (cosh(-rapidity2) + q * cosh(rapidity)); // x_B - x_A
    double argument_yB_yA = -impact_parameter;   // y_B - y_A
    double profile_12 = weight.profile_chi(argument_xB_xA, argument_yB_yA, 0.,
                                           radius_width1); // w_A(x_B)

    double argument_xA_xB =
        (separation / (q + 1)) *
        (-q * cosh(rapidity) - cosh(-rapidity2)); // x_A - x_B
    double argument_yA_yB = impact_parameter;     // y_A - y_B
    double profile_21 = weight.profile_chi(argument_xA_xB, argument_yA_yB, 0.,
                                           radius_width2); // w_B(x_A)
    double profile_22 =
        weight.profile_chi(0., 0., 0., radius_width2); // w_B(x_B)

    // h_A and h_B from Eq.(22)
    double value1 = (-profile_21 * delta_2 + profile_22 * delta_1) /
                    (profile_11 * profile_22 - profile_12 * profile_21);
    double value2 = (profile_11 * delta_2 - profile_12 * delta_1) /
                    (profile_11 * profile_22 - profile_12 * profile_21);

    // Correct the conformal factor
    double chi_corrected = chi_plain + profile1 * value1 + profile2 * value2;
    vars.chi = pow(chi_corrected, -4.0 / conformal_power);

    // Define initial lapse
    vars.lapse += sqrt(star1_vars.lapse * star1_vars.lapse +
                       star2_vars.lapse * star2_vars.lapse - 1.);

    // Define initial trace of K and A_ij
    double one_third = 1. / 3.;
    FOR2(i, j)
    vars.h[i][j] = pow(chi_plain, -4.0 / conformal_power) * gammaLL[i][j];
    FOR4(i, j, k, l)
    KLL[i][j] += gammaLL[i][l] * (gammaUU_1[l][k] * star1_vars.KLL[k][j] +
                                  gammaUU_2[l][k] * star2_vars.KLL[k][j]);
    FOR2(i, j) vars.K += KLL[i][j] * gammaUU[i][j];
    FOR2(i, j)
    vars.A[i][j] = pow(chi_plain, -4.0 / conformal_power) *
                   (KLL[i][j] - one_third * vars.K * gammaLL[i][j]);

    current_cell.store_vars(vars);
}

#endif /* BINARYUNEQUALMASSFIX_IMPL_HPP_ */
