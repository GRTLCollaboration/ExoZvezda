/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BINARYPLAINSUPERPOSITION_HPP_)
#error                                                                         \
    "This file should only be included through BinaryPlainSuperposition.impl.hpp"
#endif

#ifndef BINARYPLAINSUPERPOSITION_IMPL_HPP_
#define BINARYPLAINSUPERPOSITION_IMPL_HPP_

#include "BosonStarSolution.hpp"

inline BinaryPlainSuperposition::BinaryPlainSuperposition(
    BosonStar_params_t a_params_BosonStar,
    BosonStar_params_t a_params_BosonStar2,
    Potential::params_t a_params_potential, double a_dx)
    : m_dx(a_dx), m_params_BosonStar(a_params_BosonStar),
      m_params_BosonStar2(a_params_BosonStar2),
      m_params_potential(a_params_potential)
{
}

void BinaryPlainSuperposition::compute_1d_solution(const double max_r)
{
    try
    {
        // set initial parameters and then run the solver (didnt put it in the
        // constructor)

        pout() << "Setting initial conditions for Star 1" << endl;
        m_1d_sol.set_initialcondition_params(m_params_BosonStar,
                                             m_params_potential, max_r);
        pout() << "Running the solver for Star 1" << endl;
        m_1d_sol.main();
        pout() << "Completed for star 1" << endl;

        pout() << "Setting initial conditions for Star 2" << endl;
        m_1d_sol2.set_initialcondition_params(m_params_BosonStar2,
                                              m_params_potential, max_r);
        pout() << "Running the solver for Star 2" << endl;
        m_1d_sol2.main();
        pout() << "Completed for star 2" << endl;
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class data_t>
void BinaryPlainSuperposition::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)

    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below

    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_BosonStar.star_centre);

    // Import BS parameters and option of whether this is a BS binary or BS-BH
    // binary
    double rapidity = m_params_BosonStar.BS_rapidity;
    double rapidity2 = m_params_BosonStar2.BS_rapidity;
    bool antiboson = m_params_BosonStar.antiboson;
    double separation = m_params_BosonStar.BS_separation;
    double impact_parameter = m_params_BosonStar.BS_impact_parameter;
    double q = m_params_BosonStar.mass_ratio;

    // Define boosts and coordinate objects, suppose star 1 is on the left of
    // the centre of mass and star 2 is on the right of centre of mass

    // e.g. taking the centre of mass to be the origin, then STAR2 --------
    // (origin) -------- STAR1

    // First star positioning
    double c_ = cosh(rapidity);
    double s_ = sinh(rapidity);
    double v_ = tanh(rapidity);
    double t =
        (coords.x - q * separation / (q + 1.)) * s_; // set /tilde{t} to zero
    double x = (coords.x - q * separation / (q + 1.)) * c_;
    double z = coords.z; // set /tilde{t} to zero
    double y = coords.y + q * impact_parameter / (q + 1.);
    double r = sqrt(x * x + y * y + z * z);

    // First star physical variables
    double p_ = m_1d_sol.get_A_interp(r);
    double dp_ = m_1d_sol.get_dA_interp(r);
    double omega_ = m_1d_sol.get_lapse_interp(r);
    double omega_prime_ = m_1d_sol.get_dlapse_interp(r);
    double psi_ = m_1d_sol.get_psi_interp(r);
    double psi_prime_ = m_1d_sol.get_dpsi_interp(r);

    // Get scalar field modulus, conformal factor, lapse and their gradients
    double pc_os = psi_ * psi_ * c_ * c_ - omega_ * omega_ * s_ * s_;
    double lapse_1 = omega_ * psi_ / (sqrt(pc_os));
    double lapse_2 = 1.;
    double w_ = m_1d_sol.get_w();

    // Write in phase, shift, metric componnets of star 1 and initialise metric
    // components of star 2
    double phase_ = m_params_BosonStar.phase * M_PI + w_ * t;
    double beta_x = s_ * c_ * (psi_ * psi_ - omega_ * omega_) / (pc_os);
    vars.shift[0] += beta_x;
    double g_zz_1 = psi_ * psi_;
    double g_yy_1 = psi_ * psi_;
    double g_xx_1 = pc_os;
    double g_xx_2 = 0., g_yy_2 = 0., g_zz_2 = 0., g_xx, g_yy, g_zz;

    // Add on to evolution equations
    vars.phi_Re += p_ * cos(phase_);
    vars.phi_Im += p_ * sin(phase_);
    vars.Pi_Re +=
        -(1. / lapse_1) * ((x / r) * (s_ - beta_x * c_) * dp_ * cos(phase_) -
                           w_ * (c_ - beta_x * s_) * p_ * sin(phase_));
    vars.Pi_Im +=
        -(1. / lapse_1) * ((x / r) * (s_ - beta_x * c_) * dp_ * sin(phase_) +
                           w_ * (c_ - beta_x * s_) * p_ * cos(phase_));

    // Initialise extrinsic curvature and metric with upper indices
    double KLL_1[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double KLL_2[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double KLL[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaLL[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaUU[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaUU_1[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaUU_2[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double K1 = 0., K2 = 0.;

    // Fill them in
    gammaUU_1[0][0] = 1. / g_xx_1;
    gammaUU_1[1][1] = 1. / g_yy_1;
    gammaUU_1[2][2] = 1. / g_zz_1;

    KLL_1[2][2] = -lapse_1 * s_ * x * psi_prime_ / (r * psi_);
    KLL_1[1][1] = KLL_1[2][2];
    KLL_1[0][1] = lapse_1 * c_ * s_ * (y / r) *
                  (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_1[0][2] = lapse_1 * c_ * s_ * (z / r) *
                  (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_1[1][0] = KLL_1[0][1];
    KLL_1[2][0] = KLL_1[0][2];
    KLL_1[2][1] = 0.;
    KLL_1[1][2] = 0.;
    KLL_1[0][0] = lapse_1 * (x / r) * s_ * c_ * c_ *
                  (psi_prime_ / psi_ - 2. * omega_prime_ / omega_ +
                   v_ * v_ * omega_ * omega_prime_ * pow(psi_, -2));
    FOR2(i, j) K1 += gammaUU_1[i][j] * KLL_1[i][j];

    // Here is the conformal factor that will be differently calculated
    // depending on the choice of the initial data
    double chi_;
    double chi_plain;

    // This is the effect of object 1 on object 2 and hence represents the value
    // to be substracted in the initial data from the position of object 2
    double t_p = (-separation) * s_; // set /tilde{t} to zero
    double x_p = (-separation) * c_;
    double z_p = 0.; // set /tilde{t} to zero
    double y_p = impact_parameter;
    double r_p = sqrt(x_p * x_p + y_p * y_p + z_p * z_p);
    double p_p = m_1d_sol.get_A_interp(r_p);
    double dp_p = m_1d_sol.get_dA_interp(r_p);
    double omega_p = m_1d_sol.get_lapse_interp(r_p);
    double omega_prime_p = m_1d_sol.get_dlapse_interp(r_p);
    double psi_p = m_1d_sol.get_psi_interp(r_p);
    double psi_prime_p = m_1d_sol.get_dpsi_interp(r_p);
    double pc_os_p = psi_p * psi_p * c_ * c_ - omega_p * omega_p * s_ * s_;

    // Second star positioning
    c_ = cosh(-rapidity2);
    s_ = sinh(-rapidity2);
    v_ = tanh(-rapidity2);
    t = (coords.x + separation / (q + 1.)) * s_; // set /tilde{t} to zero
    x = (coords.x + separation / (q + 1.)) * c_;
    z = coords.z;
    y = coords.y - impact_parameter / (q + 1.);
    r = sqrt(x * x + y * y + z * z);

    // Second star physical variables
    p_ = m_1d_sol2.get_A_interp(r);
    dp_ = m_1d_sol2.get_dA_interp(r);
    omega_ = m_1d_sol2.get_lapse_interp(r);
    omega_prime_ = m_1d_sol2.get_dlapse_interp(r);
    psi_ = m_1d_sol2.get_psi_interp(r);
    psi_prime_ = m_1d_sol2.get_dpsi_interp(r);

    pc_os = psi_ * psi_ * c_ * c_ - omega_ * omega_ * s_ * s_;
    lapse_2 = omega_ * psi_ / (sqrt(pc_os));
    if (antiboson)
    {
        w_ = -m_1d_sol2.get_w();
    }
    else
    {
        w_ = m_1d_sol2.get_w();
    }
    phase_ = w_ * t;
    beta_x = s_ * c_ * (psi_ * psi_ - omega_ * omega_) / (pc_os);
    vars.shift[0] += beta_x;
    g_zz_2 = psi_ * psi_;
    g_yy_2 = psi_ * psi_;
    g_xx_2 = pc_os;
    gammaUU_2[0][0] = 1. / g_xx_2;
    gammaUU_2[1][1] = 1. / g_yy_2;
    gammaUU_2[2][2] = 1. / g_zz_2;

    vars.phi_Re += p_ * cos(phase_);
    vars.phi_Im += p_ * sin(phase_);
    vars.Pi_Re +=
        -(1. / lapse_2) * ((x / r) * (s_ - beta_x * c_) * dp_ * cos(phase_) -
                           w_ * (c_ - beta_x * s_) * p_ * sin(phase_));
    vars.Pi_Im +=
        -(1. / lapse_2) * ((x / r) * (s_ - beta_x * c_) * dp_ * sin(phase_) +
                           w_ * (c_ - beta_x * s_) * p_ * cos(phase_));

    KLL_2[2][2] = -lapse_2 * s_ * x * psi_prime_ / (r * psi_);
    KLL_2[1][1] = KLL_2[2][2];
    KLL_2[0][1] = lapse_2 * c_ * s_ * (y / r) *
                  (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_2[0][2] = lapse_2 * c_ * s_ * (z / r) *
                  (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_2[1][0] = KLL_2[0][1];
    KLL_2[2][0] = KLL_2[0][2];
    KLL_2[2][1] = 0.;
    KLL_2[1][2] = 0.;
    KLL_2[0][0] = lapse_2 * (x / r) * s_ * c_ * c_ *
                  (psi_prime_ / psi_ - 2. * omega_prime_ / omega_ +
                   v_ * v_ * omega_ * omega_prime_ * pow(psi_, -2));
    FOR2(i, j) K2 += gammaUU_2[i][j] * KLL_2[i][j];

    g_xx = g_xx_1 + g_xx_2 - 1.0;
    g_yy = g_yy_1 + g_yy_2 - 1.0;
    g_zz = g_zz_1 + g_zz_2 - 1.0;

    // Now, compute upper and lower components
    gammaLL[0][0] = g_xx;
    gammaLL[1][1] = g_yy;
    gammaLL[2][2] = g_zz;
    gammaUU[0][0] = 1. / g_xx;
    gammaUU[1][1] = 1. / g_yy;
    gammaUU[2][2] = 1. / g_zz;

    // Define initial conformal factor
    chi_ = pow(g_xx * g_yy * g_zz, -1. / 3.);

    vars.chi = chi_;

    // Define initial lapse
    vars.lapse += sqrt(lapse_1 * lapse_1 + lapse_2 * lapse_2 - 1.);

    // Define initial trace of K and A_ij
    double one_third = 1. / 3.;
    FOR2(i, j) vars.h[i][j] = vars.chi * gammaLL[i][j];
    FOR4(i, j, k, l)
    KLL[i][j] += gammaLL[i][l] * (gammaUU_1[l][k] * KLL_1[k][j] +
                                  gammaUU_2[l][k] * KLL_2[k][j]);
    FOR2(i, j) vars.K += KLL[i][j] * gammaUU[i][j];
    FOR2(i, j)
    vars.A[i][j] = vars.chi * (KLL[i][j] - one_third * vars.K * gammaLL[i][j]);

    current_cell.store_vars(vars);
}

#endif /* BINARYPLAINSUPERPOSITION_IMPL_HPP_ */
