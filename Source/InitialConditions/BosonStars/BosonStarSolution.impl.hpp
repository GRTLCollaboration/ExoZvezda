/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARSOLUTION_HPP_)
#error "This file should only be included through BosonStarSolution.hpp"
#endif

#ifndef BOSONSTARSOLUTION_IMPL_HPP_
#define BOSONSTARSOLUTION_IMPL_HPP_

#include <fstream>
#include <iomanip>
#include <iostream>

BosonStarSolution::BosonStarSolution() {}

void BosonStarSolution::main()
{
    for (int iter = 0; iter < 17; iter++)
    {

        pout() << "-----------------------------------------------" << endl;
        pout() << "I am running iteration # " << iter << endl;
        pout() << "-----------------------------------------------" << endl;

        // Set the initial conditions
        initialise();

        // The bosonic frequency of a ground state almost inevitably lies
        // between 0 and 1. Check this by finding the number of zero crossings
        // for omega ansatz of unity. If we have a zero crossing, then the upper
        // bound will be omega_upper = 1.
        omega_ansatz = find_upper_omega();
        // Apply bisection algorithm to find the right frequency of the BS
        omega_true = bisect_omega(0, omega_ansatz);
        // Determine the matching index
        matching_index = find_matching_index();

        if (BS_verbosity)
        {
            pout() << "I will use the matching index of " << matching_index
                   << " and the matching radius at " << matching_index * dx
                   << endl;
        }

        // Force the scalar field to zero after the point the amplitude diverges
        force_flat(matching_index);
        rk4_asymp(matching_index, false, omega_true);
        rk4_asymp(
            matching_index, true,
            omega_true); // (true) uses large radius adaptive stepsize to get
                         // asymptotics (integrates vacuum metric out to huge
                         // radius ~ 10^10 to measure asymptotics).

        PSI_INF = psi[gridsize - 1];
        OM_INF = omega[gridsize - 1];

        if (BS_verbosity)
        {
            pout() << "PSI_INF " << PSI_INF << endl;
            pout() << "OM_INF " << OM_INF << endl;
        }

        // if (fabs(PSI_INF - 1.) < 1e-05)
        // {
        // 	double domega = par.omega0 - 1e-14;

        //   	converged = calcExtAsym(domega, par.phigrav0);

        //   	double criterion_A_v2 = calculate_criterion_A(converged,
        //   domega);

        //   	criterion_derivative_A = (criterion_A_v2 - criterion_A) /
        //   (-1e-14);

        //   	printf("criterion derivative A now: %g \n",
        //   criterion_derivative_A);

        //   	printf("The bosonic scalar field is not matched smoothly enough.
        //   Updating omega...\n");
        // }

        // Have some stopping criterion for hard-to-compute solutions.
        // Useful, when you may want to iterate through the for loop a few
        // times.
        if (fabs(PSI_INF - 1.) < 1e-05 && fabs(OM_INF - 1.) < 1e-05)
        {
            pout() << "Central Density : " << A[0] << ", PSI0 : " << psi[0]
                   << ", OM0 : " << omega[0] << ", w : " << sqrt(omega_true)
                   << endl;
            break;
        }

        // Update initial guesses
        PSC /= PSI_INF;
        OMC /= OM_INF;

        pout() << "Central Density : " << A[0] << ", PSI0 : " << PSC
               << ", OM0 : " << OMC << ", w : " << sqrt(omega_true) << endl;

        if (iter == 16)
        {
            pout() << "Ooopsies... I have reached 16 iterations now, did not "
                      "find a BS solution  ¯|_(ツ)_|¯ "
                   << endl;
            exit(0);
        }
    }

    pout() << "-----------------------------------------------" << endl;
    pout() << "Computing the diagnostics of the solution" << endl;
    pout() << "-----------------------------------------------" << endl;
    rk4_match(matching_index, false, omega_true);

    // Calculate the aspect mass and the ADM mass at the boundary of the
    // physical domain. They shoud be similar, but not equal!
    calculate_aspect_mass();
    calculate_adm_mass();
    radius = calculate_radius();
    compactness_value = boson_mass[gridsize - 1] / radius;

    pout() << "-----------------------------------------------" << endl;
    pout() << "Central Density : " << A[0] << endl;
    pout() << "ADM mass : " << adm_mass[gridsize - 1] << endl;
    pout() << "Aspect mass : " << boson_mass[gridsize - 1] << endl;
    pout() << "Radius : " << radius << endl;
    pout() << "Compactness : " << compactness_value << endl;
    pout() << "W : " << sqrt(omega_true) << endl;
}

// Initialise the 5 filed variables with their central values
void BosonStarSolution::initialise()
{
    A[0] = A0;            // initial central amplitude
    omega[0] = OMC;       // initial lapse
    psi[0] = PSC;         // initial conformal factor
    dA[0] = 0.;           // amplitude derivate
    dpsi[0] = 0.;         // conformal factor derivative
    radius_array[0] = 0.; // isotropic radius
}

// If the scalar field diverges, then we rounds down the amplitude solution to a
// value where it starts to diverge. This does not affect the zero crossings but
// allows to systematically deal with infinite and diverging values.
void BosonStarSolution::truncate_solution()
{
    bool diverges = false; // a flag to see whether the amplitude gets twice as
                           // large as the central value
    double truncation;     // value of the scalar field at which we truncate the
                           // solution

    for (int i = 0; i < gridsize; ++i)
    {
        if (!diverges)
        {
            if (fabs(A[i]) > 1.01 * A[0])
            {
                diverges = true;
                truncation = A[i];
            }
        }
        if (diverges)
        {
            A[i] = truncation;
        }
    }
}

// sets scalar field and gradient to zero after the point decided by function
// find_midint
void BosonStarSolution::force_flat(const int iter_crit)
{
    for (int i = iter_crit + 1; i < gridsize; ++i)
    {
        A[i] = 0.;
        dA[i] = 0.;
    }
}

// Find the index (radius) where the amplitude starts to diverge, this be used
// for matching in the asymptotic regime
int BosonStarSolution::find_matching_index()
{
    int matching_index;

    for (int i = 0; i < gridsize - 1; ++i)
    {
        if (fabs(A[i]) < fabs(A[i + 1]) && fabs(A[i]) < fabs(A[i - 1]))
        {
            matching_index = iofr(radius_array[i] * 0.9);
            return matching_index;
        }
    }

    return gridsize - 1;
}

// Find the largest index for a give radius r
int BosonStarSolution::iofr(double rtarget)
{
    int i;

    i = 0;
    for (i = 1; i < gridsize - 1; i++)
    {
        if (radius_array[i - 1] <= rtarget && radius_array[i] > rtarget)
        {
            break;
        }
    }

    return i;
}

// Find the upper value of omega to be used in the bisection algorithm -- for
// groudn states this is usually upper_omega = 1.
double BosonStarSolution::find_upper_omega()
{
    bool crossed;
    double omega = 1.; // ansatz
    while (true)
    {
        // initialise();
        rk4(omega); // integrate
        truncate_solution();
        crossed = found_zero_crossing();
        if (crossed)
            return omega;
        omega *= 2.;
    }
}

// Bisection algorithm for narrowing down the frequency solution.
// The bisection checks for zero crossings and narrows down the range
// [lower_omega;upper omega] until a desired difference between upper_omega and
// lower_omega is reached (user specifiable).
double BosonStarSolution::bisect_omega(double omega_min, double omega_max)
{
    int iter = 0;
    double lower_omega, middle_omega, upper_omega;
    bool crossed;

    lower_omega = omega_min;
    upper_omega = omega_max;

    while (true)
    {
        iter++;
        middle_omega = 0.5 * (lower_omega + upper_omega);
        //   initialise();
        rk4(middle_omega);
        truncate_solution();
        crossed = found_zero_crossing();
        if (crossed)
        {
            upper_omega = middle_omega; // if crossed, then the upper_omega is
                                        // too large, set it to middle_omega
        }
        else
        {
            lower_omega = middle_omega; // if did not cross, then lower_omega
                                        // may be increased to middle_omega
        }
        if (fabs(upper_omega - lower_omega) < omega_tolerance)
            return upper_omega;
        if (upper_omega == lower_omega)
            return upper_omega;
        if (iter > 100)
        {
            if (BS_verbosity)
            {
                pout()
                    << "BosonStarSolution::bisect_omega -- I have reached the "
                       "maximum number of iterations in the bisection algorithm"
                    << endl;
            }
            return upper_omega;
        }
    }
}

// Check is the amplitudes has a zero crossing
int BosonStarSolution::found_zero_crossing()
{
    int any_zero_crossings;

    for (int i = 1; i < gridsize; i++)
    {
        if (A[i] * A[i + 1] < 0.)
        {
            any_zero_crossings = 1;
            break;
        }
        if (A[i] > A[i - 1])
        {
            any_zero_crossings = 0;
            break;
        }
    }

    return any_zero_crossings;
}

// Integrate the full ODE system from r=0 to dx*gridsize.
// The solution almost inevitable blows up, but we remedy this by matching the
// asumptotics later in the main().
void BosonStarSolution::rk4(const double ww_)
{
    double k1 = 0, k2 = 0, k3 = 0, k4 = 0, q1 = 0, q2 = 0, q3 = 0, q4 = 0,
           x_ = 0., h = dx / 2.;
    const double DX = dx;
    double DX_ = DX;
    double o1, o2, o3, o4, s1, s2, s3, s4, r1, r2, r3, r4;
    int index, jmax = 0;
    radius_array[0] = 0.;

    for (int i = 1; i < gridsize; ++i)
    {
        DX_ = DX;
        jmax = 0;
        if (i < adaptive_buffer)
        {
            jmax = adaptive_stepsize_repetitions;
        }
        for (int j = 0; j <= jmax; j++)
        {
            DX_ = DX / ((double)(1 + jmax));
            h = DX_ / 2.;
            x_ = (i - 1) * dx + j * DX_;

            k1 = DX_ * A_RHS(x_, A[i - 1], dA[i - 1], psi[i - 1], dpsi[i - 1],
                             omega[i - 1], ww_);
            q1 = DX_ * DA_RHS(x_, A[i - 1], dA[i - 1], psi[i - 1], dpsi[i - 1],
                              omega[i - 1], ww_);
            o1 = DX_ * OMEGA_RHS(x_, A[i - 1], dA[i - 1], psi[i - 1],
                                 dpsi[i - 1], omega[i - 1], ww_);
            s1 = DX_ * PSI_RHS(x_, A[i - 1], dA[i - 1], psi[i - 1], dpsi[i - 1],
                               omega[i - 1], ww_);
            r1 = DX_ * DPSI_RHS(x_, A[i - 1], dA[i - 1], psi[i - 1],
                                dpsi[i - 1], omega[i - 1], ww_);

            k2 = DX_ * A_RHS(x_ + h, A[i - 1] + k1 / 2., dA[i - 1] + q1 / 2.,
                             psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                             omega[i - 1] + o1 / 2., ww_);
            q2 = DX_ * DA_RHS(x_ + h, A[i - 1] + k1 / 2., dA[i - 1] + q1 / 2.,
                              psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                              omega[i - 1] + o1 / 2., ww_);
            o2 =
                DX_ * OMEGA_RHS(x_ + h, A[i - 1] + k1 / 2., dA[i - 1] + q1 / 2.,
                                psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                                omega[i - 1] + o1 / 2., ww_);
            s2 = DX_ * PSI_RHS(x_ + h, A[i - 1] + k1 / 2., dA[i - 1] + q1 / 2.,
                               psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                               omega[i - 1] + o1 / 2., ww_);
            r2 = DX_ * DPSI_RHS(x_ + h, A[i - 1] + k1 / 2., dA[i - 1] + q1 / 2.,
                                psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                                omega[i - 1] + o1 / 2., ww_);

            k3 = DX_ * A_RHS(x_ + h, A[i - 1] + k2 / 2., dA[i - 1] + q2 / 2.,
                             psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                             omega[i - 1] + o2 / 2., ww_);
            q3 = DX_ * DA_RHS(x_ + h, A[i - 1] + k2 / 2., dA[i - 1] + q2 / 2.,
                              psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                              omega[i - 1] + o2 / 2., ww_);
            o3 =
                DX_ * OMEGA_RHS(x_ + h, A[i - 1] + k2 / 2., dA[i - 1] + q2 / 2.,
                                psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                                omega[i - 1] + o2 / 2., ww_);
            s3 = DX_ * PSI_RHS(x_ + h, A[i - 1] + k2 / 2., dA[i - 1] + q2 / 2.,
                               psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                               omega[i - 1] + o2 / 2., ww_);
            r3 = DX_ * DPSI_RHS(x_ + h, A[i - 1] + k2 / 2., dA[i - 1] + q2 / 2.,
                                psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                                omega[i - 1] + o2 / 2., ww_);

            k4 = DX_ * A_RHS(x_ + 2. * h, A[i - 1] + k3, dA[i - 1] + q3,
                             psi[i - 1] + s3, dpsi[i - 1] + r3,
                             omega[i - 1] + o3, ww_);
            q4 = DX_ * DA_RHS(x_ + 2. * h, A[i - 1] + k3, dA[i - 1] + q3,
                              psi[i - 1] + s3, dpsi[i - 1] + r3,
                              omega[i - 1] + o3, ww_);
            o4 = DX_ * OMEGA_RHS(x_ + 2. * h, A[i - 1] + k3, dA[i - 1] + q3,
                                 psi[i - 1] + s3, dpsi[i - 1] + r3,
                                 omega[i - 1] + o3, ww_);
            s4 = DX_ * PSI_RHS(x_ + 2. * h, A[i - 1] + k3, dA[i - 1] + q3,
                               psi[i - 1] + s3, dpsi[i - 1] + r3,
                               omega[i - 1] + o3, ww_);
            r4 = DX_ * DPSI_RHS(x_ + 2. * h, A[i - 1] + k3, dA[i - 1] + q3,
                                psi[i - 1] + s3, dpsi[i - 1] + r3,
                                omega[i - 1] + o3, ww_);

            index = i - 1;
            if (j == jmax)
            {
                index = i;
            }

            A[index] = A[i - 1] + (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
            dA[index] = dA[i - 1] + (q1 + 2. * q2 + 2. * q3 + q4) / 6.;
            psi[index] = psi[i - 1] + (s1 + 2. * s2 + 2. * s3 + s4) / 6.;
            dpsi[index] = dpsi[i - 1] + (r1 + 2. * r2 + 2. * r3 + r4) / 6.;
            omega[index] = omega[i - 1] + (o1 + 2. * o2 + 2. * o3 + o4) / 6.;
        }
        radius_array[i] = dx * i;
    }
}

// takes an integrated ODE system and starts at point (iter) and re-integrates
// but enforcing scalara field to decaay or be in vacuum the integral is
// adaptive in that it aaccelerates ar later radius in order to find correct
// asymptotic behaviour. It will shout if the radius reached is below 10e10 bool
// adaaptive is true if stepsize is supposed to be adaptive aat large radius and
// false for constant stepsize
void BosonStarSolution::rk4_asymp(const int iter, const bool adaptive,
                                  const double ww_)
{
    double k1 = 0, k2 = 0, k3 = 0, k4 = 0, q1 = 0, q2 = 0, q3 = 0, q4 = 0,
           x_ = iter * dx, h, delta = (double)gridsize;
    const double DX = dx;
    double DX_ = DX;
    double o1, o2, o3, o4, s1, s2, s3, s4, r1, r2, r3, r4;
    double N_ = gridsize - iter, L_ = pow(9., 9);
    int i_;

    double k_ = log(L_) / N_;

    for (int i = iter + 1; i < gridsize; ++i)
    {
        i_ = double(i - iter);
        if (adaptive)
        {
            if (x_ < 8e8)
            {
                DX_ = (exp(k_) - 1.) * exp(k_ * i_);
            }
            else
            {
                DX_ = DX;
            }
        }
        h = DX_ / 2.;

        // k1 =
        // dx*small_P_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);
        // q1 = dx*DP_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);
        o1 = DX_ * OMEGA_RHS(x_, A[i - 1], dA[i - 1], psi[i - 1], dpsi[i - 1],
                             omega[i - 1], ww_);
        s1 = DX_ * PSI_RHS(x_, A[i - 1], dA[i - 1], psi[i - 1], dpsi[i - 1],
                           omega[i - 1], ww_);
        r1 = DX_ * DPSI_RHS(x_, A[i - 1], dA[i - 1], psi[i - 1], dpsi[i - 1],
                            omega[i - 1], ww_);

        // k2 = dx*small_P_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] +
        // s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_); q2 = dx*DP_RHS(x_ +
        // h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] +
        // r1/2.,omega[i-1] + o1/2.,ww_);
        o2 = DX_ * OMEGA_RHS(x_ + h, A[i - 1] + k1 / 2., dA[i - 1] + q1 / 2.,
                             psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                             omega[i - 1] + o1 / 2., ww_);
        s2 = DX_ * PSI_RHS(x_ + h, A[i - 1] + k1 / 2., dA[i - 1] + q1 / 2.,
                           psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                           omega[i - 1] + o1 / 2., ww_);
        r2 = DX_ * DPSI_RHS(x_ + h, A[i - 1] + k1 / 2., dA[i - 1] + q1 / 2.,
                            psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                            omega[i - 1] + o1 / 2., ww_);

        // k3 = dx*small_P_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] +
        // s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_); q3 = dx*DP_RHS(x_ +
        // h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] +
        // r2/2.,omega[i-1] + o2/2.,ww_);
        o3 = DX_ * OMEGA_RHS(x_ + h, A[i - 1] + k2 / 2., dA[i - 1] + q2 / 2.,
                             psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                             omega[i - 1] + o2 / 2., ww_);
        s3 = DX_ * PSI_RHS(x_ + h, A[i - 1] + k2 / 2., dA[i - 1] + q2 / 2.,
                           psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                           omega[i - 1] + o2 / 2., ww_);
        r3 = DX_ * DPSI_RHS(x_ + h, A[i - 1] + k2 / 2., dA[i - 1] + q2 / 2.,
                            psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                            omega[i - 1] + o2 / 2., ww_);

        // k4 = dx*small_P_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] +
        // s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_); q4 = dx*DP_RHS(x_
        // + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] +
        // r3,omega[i-1] + o3,ww_);
        o4 = DX_ * OMEGA_RHS(x_ + 2. * h, A[i - 1] + k3, dA[i - 1] + q3,
                             psi[i - 1] + s3, dpsi[i - 1] + r3,
                             omega[i - 1] + o3, ww_);
        s4 = DX_ * PSI_RHS(x_ + 2. * h, A[i - 1] + k3, dA[i - 1] + q3,
                           psi[i - 1] + s3, dpsi[i - 1] + r3, omega[i - 1] + o3,
                           ww_);
        r4 = DX_ * DPSI_RHS(x_ + 2. * h, A[i - 1] + k3, dA[i - 1] + q3,
                            psi[i - 1] + s3, dpsi[i - 1] + r3,
                            omega[i - 1] + o3, ww_);

        A[i] = 0.;  // p[i-1] + (k1 + 2.*k2 + 2.*k3 + k4)/6.;
        dA[i] = 0.; // dp[i-1] + (q1 + 2.*q2 + 2.*q3 + q4)/6.;
        psi[i] = psi[i - 1] + (s1 + 2. * s2 + 2. * s3 + s4) / 6.;
        dpsi[i] = dpsi[i - 1] + (r1 + 2. * r2 + 2. * r3 + r4) / 6.;
        omega[i] = omega[i - 1] + (o1 + 2. * o2 + 2. * o3 + o4) / 6.;
        x_ += DX_;
        if (!adaptive)
        {
            radius_array[i] = i * dx;
        }
    }

    if (adaptive and x_ < 8e7)
    {
        pout() << "\33[30;41m"
               << " Asymptotic Radius Too Small"
               << "\x1B[0m" << endl;
        pout() << x_ << endl;
    }
}

// takes an integrated ODE system and starts at point (iter) and re-integrates
// but enforcing scalara field to decaay or be in vacuum the integral is
// adaptive in that it aaccelerates ar later radius in order to find correct
// asymptotic behaviour. It will shout if the radius reached is below 10e10 bool
// adaaptive is true if stepsize is supposed to be adaptive aat large radius and
// false for constant stepsize
void BosonStarSolution::rk4_match(const int iter, const bool adaptive,
                                  const double ww_)
{
    double x_, h;
    const double DX = dx;
    double DX_ = DX;
    double o1, o2, o3, o4, s1, s2, s3, s4, r1, r2, r3, r4;
    double N_ = gridsize - iter, L_ = pow(9., 9);
    int i_;

    double k_ = log(L_) / N_;

    double r, dr;
    double c1, c2;
    double Amp, eta, mass, arealr, om0, epsilon;

    om0 = pow(1 / omega[iter], 2);
    mass = -psi[iter] * dpsi[iter] * radius_array[iter] * radius_array[iter];
    epsilon = mass * (1 - 2 * ww_) / sqrt(1 - ww_ / om0);

    r = radius_array[iter];
    arealr = r + mass + mass * mass / (4 * r);

    c1 = A[iter] * pow(r, 1 + epsilon) * exp(r * sqrt(1 - ww_ / om0));
    c2 = dA[iter] * pow(r, 1 + epsilon) * exp(r * sqrt(1 - ww_ / om0));

    if (BS_verbosity)
    {
        pout() << "Constant A " << c1 << endl;
        pout() << "Constant B " << c2 << endl;
    }

    for (int i = iter + 1; i < gridsize; ++i)
    {
        dr = radius_array[i] - radius_array[i - 1];

        h = DX_ / 2.;

        // 1st RK step
        r = radius_array[i - 1];
        om0 = pow(1 / omega[i - 1], 2);
        mass = -psi[i - 1] * dpsi[i - 1] * radius_array[i - 1] *
               radius_array[i - 1];
        epsilon = mass * (1 - 2 * ww_) / sqrt(1 - ww_ / om0);
        arealr = r + mass + mass * mass / (4 * r);
        Amp = c1 * exp(-r * sqrt(1 - ww_ / om0)) * pow(r, -1 - epsilon);
        eta = c2 * exp(-r * sqrt(1 - ww_ / om0)) * pow(r, -1 - epsilon);
        o1 = dr *
             OMEGA_RHS(r, Amp, eta, psi[i - 1], dpsi[i - 1], omega[i - 1], ww_);
        s1 = dr *
             PSI_RHS(r, Amp, eta, psi[i - 1], dpsi[i - 1], omega[i - 1], ww_);
        r1 = dr *
             DPSI_RHS(r, Amp, eta, psi[i - 1], dpsi[i - 1], omega[i - 1], ww_);

        // 2nd RK step
        r = radius_array[i - 1] + 0.5 * dr;
        om0 = pow(1 / (omega[i - 1] + o1 / 2.), 2);
        mass = -(psi[i - 1] + s1 / 2.) * (dpsi[i - 1] + r1 / 2.) * r * r;
        epsilon = mass * (1 - 2 * ww_) / sqrt(1 - ww_ / om0);
        arealr = r + mass + mass * mass / (4 * r);
        Amp = c1 * exp(-r * sqrt(1 - ww_ / om0)) * pow(r, -1 - epsilon);
        eta = c2 * exp(-r * sqrt(1 - ww_ / om0)) * pow(r, -1 - epsilon);
        o2 = dr * OMEGA_RHS(r, Amp, eta, psi[i - 1] + s1 / 2.,
                            dpsi[i - 1] + r1 / 2., omega[i - 1] + o1 / 2., ww_);
        s2 = dr * PSI_RHS(r, Amp, eta, psi[i - 1] + s1 / 2.,
                          dpsi[i - 1] + r1 / 2., omega[i - 1] + o1 / 2., ww_);
        r2 = dr * DPSI_RHS(r, Amp, eta, psi[i - 1] + s1 / 2.,
                           dpsi[i - 1] + r1 / 2., omega[i - 1] + o1 / 2., ww_);

        // 3rd RK step
        r = radius_array[i - 1] + 0.5 * dr;
        om0 = pow(1 / (omega[i - 1] + o2 / 2.), 2);
        mass = -(psi[i - 1] + s2 / 2.) * (dpsi[i - 1] + r2 / 2.) * r * r;
        epsilon = mass * (1 - 2 * ww_) / sqrt(1 - ww_ / om0);
        arealr = r + mass + mass * mass / (4 * r);
        Amp = c1 * exp(-r * sqrt(1 - ww_ / om0)) * pow(r, -1 - epsilon);
        eta = c2 * exp(-r * sqrt(1 - ww_ / om0)) * pow(r, -1 - epsilon);
        o3 = dr * OMEGA_RHS(r, Amp, eta, psi[i - 1] + s2 / 2.,
                            dpsi[i - 1] + r2 / 2., omega[i - 1] + o2 / 2., ww_);
        s3 = dr * PSI_RHS(r, Amp, eta, psi[i - 1] + s2 / 2.,
                          dpsi[i - 1] + r2 / 2., omega[i - 1] + o2 / 2., ww_);
        r3 = dr * DPSI_RHS(r, Amp, eta, psi[i - 1] + s2 / 2.,
                           dpsi[i - 1] + r2 / 2., omega[i - 1] + o2 / 2., ww_);

        // 4th RK step
        r = radius_array[i];
        om0 = pow(1 / (omega[i - 1] + o3), 2);
        mass = -(psi[i - 1] + s3) * (dpsi[i - 1] + r3) * r * r;
        epsilon = mass * (1 - 2 * ww_) / sqrt(1 - ww_ / om0);
        arealr = r + mass + mass * mass / (4 * r);
        Amp = c1 * exp(-r * sqrt(1 - ww_ / om0)) * pow(r, -1 - epsilon);
        eta = c2 * exp(-r * sqrt(1 - ww_ / om0)) * pow(r, -1 - epsilon);
        o4 = dr * OMEGA_RHS(r, Amp, eta, psi[i - 1] + s3, dpsi[i - 1] + r3,
                            omega[i - 1] + o3, ww_);
        s4 = dr * PSI_RHS(r, Amp, eta, psi[i - 1] + s3, dpsi[i - 1] + r3,
                          omega[i - 1] + o3, ww_);
        r4 = dr * DPSI_RHS(r, Amp, eta, psi[i - 1] + s3, dpsi[i - 1] + r3,
                           omega[i - 1] + o3, ww_);

        // Update variables
        psi[i] = psi[i - 1] + (s1 + 2. * s2 + 2. * s3 + s4) / 6.;
        dpsi[i] = dpsi[i - 1] + (r1 + 2. * r2 + 2. * r3 + r4) / 6.;
        omega[i] = omega[i - 1] + (o1 + 2. * o2 + 2. * o3 + o4) / 6.;
        r = radius_array[i];
        om0 = pow(1 / omega[i], 2);
        mass = -psi[i] * dpsi[i] * r * r;
        epsilon = mass * (1 - 2 * ww_) / sqrt(1 - ww_ / om0);
        arealr = r + mass + mass * mass / (4 * r);
        A[i] = c1 * exp(-r * sqrt(1 - ww_ / om0)) * pow(r, -1 - epsilon);
        dA[i] = c2 * exp(-r * sqrt(1 - ww_ / om0)) * pow(r, -1 - epsilon);
    }
}

// these functions return the right hand side of the ode's
// small_P_RHS is valid for large redius when the scalar field is small.

double BosonStarSolution::small_A_RHS(const double x, const double A,
                                      const double DA, const double PSI,
                                      const double DPSI, const double OM,
                                      const double ww_)
{
    double RHS = -A * PSI * sqrt(DV(A) - ww_ / (OM * OM));
    return RHS;
}

double BosonStarSolution::A_RHS(const double x, const double A, const double DA,
                                const double PSI, const double DPSI,
                                const double OM, const double ww_)
{
    double RHS = DA;
    return RHS;
}

double BosonStarSolution::DA_RHS(const double x, const double A,
                                 const double DA, const double PSI,
                                 const double DPSI, const double OM,
                                 const double ww_)
{
    double r = ((x == 0.) ? eps : x);
    double DOM = OMEGA_RHS(x, A, DA, PSI, DPSI, OM, ww_);
    return A * PSI * PSI * (DV(A) - ww_ / (OM * OM)) -
           DA * (DOM / OM + DPSI / PSI + 2. / r);
}

double BosonStarSolution::PSI_RHS(const double x, const double A,
                                  const double DA, const double PSI,
                                  const double DPSI, const double OM,
                                  const double ww_)
{
    double RHS = DPSI;
    return RHS;
}

double BosonStarSolution::DPSI_RHS(const double x, const double A,
                                   const double DA, const double PSI,
                                   const double DPSI, const double OM,
                                   const double ww_)
{
    double r = ((x == 0.) ? eps : x);
    return 0.5 * DPSI * DPSI / PSI - 2. * DPSI / r -
           2. * M_PI * PSI *
               (PSI * PSI * V(A) + DA * DA +
                ww_ * A * A * PSI * PSI / (OM * OM));
}

double BosonStarSolution::OMEGA_RHS(const double x, const double A,
                                    const double DA, const double PSI,
                                    const double DPSI, const double OM,
                                    const double ww_)
{
    double r = ((x == 0.) ? eps : x);
    return (OM / (x * DPSI + PSI)) *
           (2. * M_PI * x * PSI *
                (DA * DA - PSI * PSI * V(A) +
                 ww_ * A * A * PSI * PSI / (OM * OM)) -
            DPSI - 0.5 * x * DPSI * DPSI / PSI);
}

// V is klein gordon potential and DV is its gradient. Depends on #define
// star_type at top
double BosonStarSolution::V(const double A)
{
    if (!solitonic)
    {
        return MM * A * A + 0.5 * lambda * A * A * A * A;
    }
    else
    {
        return MM * A * A * pow((1. - 2. * pow(A / sigma, 2)), 2);
    }
}
double BosonStarSolution::DV(const double A)
{
    if (!solitonic)
    {
        return MM + lambda * A * A;
    }
    else
    {
        return MM - 8. * MM * pow(A / sigma, 2) + 12. * MM * pow(A / sigma, 4);
    }
}

void BosonStarSolution::calculate_aspect_mass()
{
    for (int i = 0; i < gridsize; ++i)
    {
        boson_mass[i] = 2. * radius_array[i] * (sqrt(psi[i]) - 1.);
    }
}

void BosonStarSolution::calculate_adm_mass()
{
    for (int i = 0; i < gridsize; ++i)
    {
        adm_mass[i] = -psi[i] * dpsi[i] * radius_array[i] * radius_array[i];
    }
}

double BosonStarSolution::calculate_radius()
{
    int i;

    for (i = gridsize - 2; i >= 0; i--)
        if (boson_mass[i] < 99.9 / 100.0 * boson_mass[gridsize - 1])
            break;

    pout() << "radius = " << radius_array[i + 1] << endl;

    return radius_array[i + 1];
}

// counts how many times the function crosses the axis
int BosonStarSolution::crossings()
{
    int number = 0;
    for (int i = 0; i < gridsize - 1; ++i)
    {
        if (A[i] * A[i + 1] < 0)
        {
            number += 1;
        }
    }
    return number;
}

// 4th order error (cubic interpolation) for field. shouts if asked to fetch a
// value outside the ode solution
double BosonStarSolution::get_A_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 =
        ((iter == 0)
             ? A[1]
             : A[iter - 1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = A[iter];
    f3 = A[iter + 1];
    f4 = A[iter + 2];

    if (iter > gridsize - 3)
    {
        pout() << "FArrayBox domain exceeding star radius!" << endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}

double BosonStarSolution::get_dA_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? dA[1] : dA[iter - 1]); // conditionl/ternary imposing
                                               // zero gradeint at r=0
    f2 = dA[iter];
    f3 = dA[iter + 1];
    f4 = dA[iter + 2];

    if (iter > gridsize - 3)
    {
        pout() << "FArrayBox domain exceeding star radius!" << endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}

double BosonStarSolution::get_lapse_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? omega[1]
                      : omega[iter - 1]); // conditionl/ternary imposing zero
                                          // gradeint at r=0
    f2 = omega[iter];
    f3 = omega[iter + 1];
    f4 = omega[iter + 2];

    if (iter > gridsize - 3)
    {
        pout() << "FArrayBox domain exceeding star radius!" << endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}

double BosonStarSolution::get_psi_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? psi[1] : psi[iter - 1]); // conditionl/ternary imposing
                                                 // zero gradeint at r=0
    f2 = psi[iter];
    f3 = psi[iter + 1];
    f4 = psi[iter + 2];

    if (iter > gridsize - 3)
    {
        pout() << "FArrayBox domain exceeding star radius!" << endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}

double BosonStarSolution::get_dpsi_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 =
        ((iter == 0) ? dpsi[1] : dpsi[iter - 1]); // conditionl/ternary imposing
                                                  // zero gradeint at r=0
    f2 = dpsi[iter];
    f3 = dpsi[iter + 1];
    f4 = dpsi[iter + 2];

    if (iter > gridsize - 3)
    {
        pout() << "FArrayBox domain exceeding star radius!" << endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}

double BosonStarSolution::get_dlapse_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? omega[1] : omega[iter - 1]);
    f2 = omega[iter];
    f3 = omega[iter + 1];
    f4 = omega[iter + 2];

    if (iter > gridsize - 3)
    {
        pout() << "FArrayBox domain exceeding star radius!" << endl;
    }

    // do the cubic spline (for gradient now), from mathematica script written
    // by Robin (rc634@cam.ac.uk)
    interpolated_value =
        (1. / (24. * dx)) *
        ((f1 - 27. * f2 + 27. * f3 - f4) + 12. * a * (f1 - f2 - f3 + f4) -
         12. * a * a * (f1 - 3. * f2 + 3. * f3 - f4));
    return interpolated_value;
}

double BosonStarSolution::get_mass() const { return boson_mass[gridsize - 1]; }

double BosonStarSolution::get_w() const { return sqrt(omega_true); }

void BosonStarSolution::set_initialcondition_params(
    BosonStar_params_t m_params_BosonStar,
    Potential::params_t m_params_potential, const double max_r)
{
    gridsize = m_params_BosonStar.gridpoints;
    adaptive_buffer =
        0.; // gridsize/10; // numer of gridpoints to intergate more carefully
    A.resize(gridsize);                  // scalar field modulus
    dA.resize(gridsize);                 // scalar field modulus gradient
    psi.resize(gridsize);                // conformal factor
    dpsi.resize(gridsize);               // conformal factor gradient
    omega.resize(gridsize);              // lapse
    radius_array.resize(gridsize);       // isotropic radius
    areal_radius_array.resize(gridsize); // arealradius
    boson_mass.resize(gridsize);         // mass of the BS
    adm_mass.resize(gridsize);           // ADM BS mass
    compactness.resize(gridsize);        // compactness of a BS

    A0 = m_params_BosonStar.central_amplitude_CSF;
    MM = m_params_potential.scalar_mass * m_params_potential.scalar_mass;
    lambda = m_params_potential.phi4_coeff;
    solitonic = m_params_potential.solitonic;
    sigma = m_params_potential.sigma_solitonic;
    BS_verbosity = m_params_BosonStar.BS_solver_verbosity;
    PSC = m_params_BosonStar.PSC;
    OMC = m_params_BosonStar.OMC;
    L = max_r * 1.05; // just to make sure the function domain is slightly
                      // larger than the required cube
    dx = L / double((gridsize - 1));
}

void BosonStarSolution::shout() const { pout() << "Haliboombah!" << endl; }

void BosonStarSolution::output_csv()
{
    std::ofstream A_file, dA_file, psi_file, dpsi_file, omega_file, r_file,
        mass_file;
    A_file.open("A.csv");
    dA_file.open("dA.csv");
    psi_file.open("psi.csv");
    dpsi_file.open("dpsi.csv");
    omega_file.open("omega.csv");
    r_file.open("r.csv");
    mass_file.open("mass.csv");

    for (int i = 0; i < gridsize; i++)
    {
        A_file << A[i] << "," << endl;
        dA_file << dA[i] << "," << endl;
        psi_file << psi[i] << "," << endl;
        dpsi_file << dpsi[i] << "," << endl;
        omega_file << omega[i] << "," << endl;
        r_file << radius_array[i] << "," << endl;
        mass_file << boson_mass[i] << "," << endl;
    }

    A_file.close();
    dA_file.close();
    psi_file.close();
    dpsi_file.close();
    omega_file.close();
    r_file.close();
    mass_file.close();
}

#endif /* BOSONSTARSOLUTION_IMPL_HPP_ */
