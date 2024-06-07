/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARSOLUTION_HPP_)
#error "This file should only be included through BosonStarSolution.hpp"
#endif

#ifndef BOSONSTARSOLUTION_IMPL_HPP_
#define BOSONSTARSOLUTION_IMPL_HPP_

BosonStarSolution::BosonStarSolution() {}

void BosonStarSolution::main()
{
    // finds the eigenvalue corresponding to intitial metric central values.
    // WW is big then it descends to the correct w+ and w- then uses interval
    // bisection to find best ww to machine precision

    WW = find_WW_soliton();
    ww = ww_IB_soliton(0, WW);
    mid_int = find_midint();

    // force the scalar field to zero after the turning point and reintegrate
    // the lapse and shift
    force_flat(mid_int);
    rk4_asymp(mid_int - 1, false, ww);
    //fix();
    rk4_asymp(mid_int, true,
              ww); // (true) uses large radius adaptive stepsize to get
                   // asymptotics (integrates vacuum metric out to huge radius ~
                   // 10^10 to measure asymptotics).
    PSI_INF = psi[gridsize - 1];
    OM_INF = omega[gridsize - 1];

    for (int q = 0; q < 2; q++)
    {
        PSC /= PSI_INF;
        OMC /= OM_INF;
        WW = find_WW_soliton();
        ww = ww_IB_soliton(0, WW);
        mid_int = find_midint();
        fix();
        force_flat(mid_int);
        rk4_asymp(mid_int - 1, false, ww);
        // calculate the ADM mass and aspect mass at the edge of physical
        // domain. Aspect is more accurate but ADM should be relativiely
        // similar.
        adm_mass = -psi[gridsize - 1] * dpsi[gridsize - 1] *
                   radius_array[gridsize - 1] * radius_array[gridsize - 1];
        aspect_mass =
            2. * radius_array[gridsize - 1] * (sqrt(psi[gridsize - 1]) - 1.);
        // std::cout << "Central Density : " << p[0] << ", ADM mass : " <<
        // adm_mass << ", aspect mass : "  << aspect_mass << ", w : "  <<
        // sqrt(ww) << std::endl;
        initialise();
        rk4_asymp(mid_int, true, ww);
        PSI_INF = psi[gridsize - 1];
        OM_INF = omega[gridsize - 1];

        /*PSC/=PSI_INF;
        OMC/=OM_INF;
        ww/=OM_INF*OM_INF;
        initialise();
        rk4(ww);
        mid_int = find_midint();
        rk4_asymp(mid_int,true,ww);
        PSI_INF = psi[gridsize-1];
        OM_INF = omega[gridsize-1];*/
    }

    rk4_asymp(mid_int - 1, false, ww);
    // calculate the ADM mass and aspect mass at the edge of physical domain.
    // Aspect is more accurate but ADM should be relativiely similar.
    adm_mass = -psi[gridsize - 1] * dpsi[gridsize - 1] *
               radius_array[gridsize - 1] * radius_array[gridsize - 1];
    aspect_mass =
        2. * radius_array[gridsize - 1] * (sqrt(psi[gridsize - 1]) - 1.);
    std::cout << "Central Density : " << p[0] << ", ADM mass : " << adm_mass
              << ", aspect mass : " << aspect_mass << ", w : " << sqrt(ww)
              << ", r99 : " << get_r(0.99) << ", mass error : "
              << fabs(2. * (aspect_mass - adm_mass) / (adm_mass + aspect_mass))
              << std::endl;
}

// initiaalises the 5 filed variables with their central values
void BosonStarSolution::initialise()
{
    p[0] = PC;
    omega[0] = OMC;
    psi[0] = PSC;
    // field gradients must be zero at zero radius for physical solutions
    dp[0] = 0.;
    dpsi[0] = 0.;
    radius_array[0] = 0.;
}

// if the scalar field diverges it rounds down to a value wiht the same sign.
// This is to not affect axis crossings function (crossings) and let it
// accurately deal with inf/nan
void BosonStarSolution::fix()
{
    bool borked = false; // turns true if function gets over twice as large as
                         // central (r=0) value
    double truncation;
    for (int i = 0; i < gridsize; ++i)
    {
        if (not borked)
        {
            if (fabs(p[i]) > 1.01 * p[0])
            {
                borked = true;
                truncation = p[i];
            }
        }
        if (borked)
        {
            p[i] = truncation;
        }
    }
}

// sets scalar field and gradient to zero after the point decided by function
// find_midint
void BosonStarSolution::force_flat(const int iter_crit)
{
    for (int i = iter_crit + 1; i < gridsize; ++i)
    {
        p[i] = 0.;
        dp[i] = 0.;
    }
}

// finds the integer index at which the scalar field needs to be truncated
int BosonStarSolution::find_midint()
{
    int crossings = 0, mid_int;

    // follow the correct amount of crossings
    for (int i = 0; i < gridsize - 1; ++i)
    {
        if (crossings == EIGEN)
        {
            mid_int = i;
            break;
        }
        if (p[i] * p[i + 1] < 0.)
        {
            crossings += 1;
        }
    }

    // climb the hill if crossings != 0, this will exit loop immediately if
    // there is no crossings
    for (int i = mid_int + 5; i < gridsize - 1; ++i)
    {
        if (fabs(p[i + 1]) < fabs(p[i]))
        {
            mid_int = i;
            break;
        }
    }

    for (int i = mid_int + 5; i < gridsize - 1; ++i)
    {
        if (fabs(p[i + 1]) > fabs(p[i]))
        {
            mid_int = i;
            // std::cout << "Truncation error: " << p[i]/PC << std::endl;
            return mid_int;
        }
    }
    return gridsize - 1;
}

// finds an eigenvalue with a lot of nodes, (20 + desired eigenstate) by defualt
double BosonStarSolution::find_WW()
{
    int eigenstate;
    double WW_ = 1.;
    while (true)
    {
        initialise();
        rk4(WW_);
        fix();
        eigenstate = crossings();
        if (eigenstate >= 20 + EIGEN)
        {
            return WW_;
        }
        WW_ *= 2.;
    }
}

double BosonStarSolution::find_WW_soliton()
{
    bool crossed;
    double WW_ = 1.;
    while (true)
    {
        initialise();
        rk4(WW_);
        fix();
        crossed = soliton_eigen();
        if (crossed)
            return WW_;
        if (WW_>10e10)
        {
            return WW_;
        }
        WW_ *= 2.;
    }
}

// calculates the lower limit for eigenvalue to be used in interval bisection
double BosonStarSolution::ww_min(const double WW_)
{
    int accuracy = 400, eigenstate;
    double ww_, lower_ww_;
    for (int i = 0; i < accuracy; i++)
    {
        ww_ = WW_ * (double)(accuracy - i) / (double)accuracy;
        initialise();
        rk4(ww_);
        fix();
        eigenstate = crossings();
        if (eigenstate <= EIGEN)
        {
            lower_ww_ = ww_;
            return lower_ww_;
        }
    }
    return 0.;
}

// calculates the upper value for eigenvalue to be used in interval bisection
double BosonStarSolution::ww_max(const double WW_, const double lower_ww_)
{
    int accuracy = 400, eigenstate;
    double ww_, upper_ww_;
    for (int i = 0; i < accuracy; i++)
    {
        ww_ = lower_ww + WW_ * ((double)i) / ((double)accuracy);
        initialise();
        rk4(ww_);
        fix();
        eigenstate = crossings();
        if (eigenstate > EIGEN)
        {
            upper_ww_ = ww_;
            return upper_ww_;
        }
    }
    return 0.;
}

// takes in an upper and lower eigenvalue and uses interval bisection to find
// the solution inbetween
double BosonStarSolution::ww_IB(double lower_ww_, double upper_ww_)
{
    int iter = 0, itermax, eigenstate, decimal_places_of_omega = 25.;
    double middle_ww_;

    itermax = (int)((log(upper_ww_) + decimal_places_of_omega * log(10.)) /
                    log(2.)); // calculate number of bisections needed (simple
                              // pen and paper calculation)
    while (true)
    {
        iter++;
        middle_ww_ = 0.5 * (upper_ww_ + lower_ww_);
        initialise();
        rk4(middle_ww_);
        fix();
        eigenstate = crossings();
        if (eigenstate > EIGEN)
        {
            upper_ww_ = middle_ww_;
        }
        else
        {
            lower_ww_ = middle_ww_;
        }
        if ((upper_ww_ - lower_ww_) < ww_tolerance)
            return upper_ww_;
        if (upper_ww_ == lower_ww_)
            return upper_ww_;
        if (iter > 65)
            return upper_ww_;
    }
}

double BosonStarSolution::ww_IB_soliton(double lower_ww_, double upper_ww_)
{
    int iter = 0, itermax, eigenstate, decimal_places_of_omega = 25.;
    double middle_ww_;
    bool crossed;

    itermax = (int)((log(upper_ww_) + decimal_places_of_omega * log(10.)) /
                    log(2.)); // calculate number of bisections needed (simple
                              // pen and paper calculation)
    while (true)
    {
        iter++;
        middle_ww_ = 0.5 * (upper_ww_ + lower_ww_);
        initialise();
        rk4(middle_ww_);
        fix();
        crossed = soliton_eigen();
        if (crossed)
        {
            upper_ww_ = middle_ww_;
        }
        else
        {
            lower_ww_ = middle_ww_;
        }
        if ((upper_ww_ - lower_ww_) < ww_tolerance)
            return upper_ww_;
        if (upper_ww_ == lower_ww_)
            return upper_ww_;
        if (iter > 100)
            return upper_ww_;
    }
}

bool BosonStarSolution::soliton_eigen()
{
    for (int i = 2; i < gridsize-6; i++)
    {
        if (p[i] * p[i + 1] < 0.)
            return true;
        if (p[i-2] * p[i + 2] < 0.)
            return true;
        if (p[i] > p[i - 1])
            return false;
    }
    return false;
}

// integrate the full ODE system from r=0 to dx*gridsize, this may (probably
// will) blow up at laarge raadius, but it is fixed by other functions
// aafterwards. it has smaller stepsize for the first (adaptive_buffer) steps.
void BosonStarSolution::rk4(const double ww_)
{
    double k1, k2, k3, k4, q1, q2, q3, q4, x_ = 0., h = dx / 2.;
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

            k1 = DX_ * P_RHS(x_, p[i - 1], dp[i - 1], psi[i - 1], dpsi[i - 1],
                             omega[i - 1], ww_);
            q1 = DX_ * DP_RHS(x_, p[i - 1], dp[i - 1], psi[i - 1], dpsi[i - 1],
                              omega[i - 1], ww_);
            o1 = DX_ * OMEGA_RHS(x_, p[i - 1], dp[i - 1], psi[i - 1],
                                 dpsi[i - 1], omega[i - 1], ww_);
            s1 = DX_ * PSI_RHS(x_, p[i - 1], dp[i - 1], psi[i - 1], dpsi[i - 1],
                               omega[i - 1], ww_);
            r1 = DX_ * DPSI_RHS(x_, p[i - 1], dp[i - 1], psi[i - 1],
                                dpsi[i - 1], omega[i - 1], ww_);

            k2 = DX_ * P_RHS(x_ + h, p[i - 1] + k1 / 2., dp[i - 1] + q1 / 2.,
                             psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                             omega[i - 1] + o1 / 2., ww_);
            q2 = DX_ * DP_RHS(x_ + h, p[i - 1] + k1 / 2., dp[i - 1] + q1 / 2.,
                              psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                              omega[i - 1] + o1 / 2., ww_);
            o2 =
                DX_ * OMEGA_RHS(x_ + h, p[i - 1] + k1 / 2., dp[i - 1] + q1 / 2.,
                                psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                                omega[i - 1] + o1 / 2., ww_);
            s2 = DX_ * PSI_RHS(x_ + h, p[i - 1] + k1 / 2., dp[i - 1] + q1 / 2.,
                               psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                               omega[i - 1] + o1 / 2., ww_);
            r2 = DX_ * DPSI_RHS(x_ + h, p[i - 1] + k1 / 2., dp[i - 1] + q1 / 2.,
                                psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                                omega[i - 1] + o1 / 2., ww_);

            k3 = DX_ * P_RHS(x_ + h, p[i - 1] + k2 / 2., dp[i - 1] + q2 / 2.,
                             psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                             omega[i - 1] + o2 / 2., ww_);
            q3 = DX_ * DP_RHS(x_ + h, p[i - 1] + k2 / 2., dp[i - 1] + q2 / 2.,
                              psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                              omega[i - 1] + o2 / 2., ww_);
            o3 =
                DX_ * OMEGA_RHS(x_ + h, p[i - 1] + k2 / 2., dp[i - 1] + q2 / 2.,
                                psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                                omega[i - 1] + o2 / 2., ww_);
            s3 = DX_ * PSI_RHS(x_ + h, p[i - 1] + k2 / 2., dp[i - 1] + q2 / 2.,
                               psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                               omega[i - 1] + o2 / 2., ww_);
            r3 = DX_ * DPSI_RHS(x_ + h, p[i - 1] + k2 / 2., dp[i - 1] + q2 / 2.,
                                psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                                omega[i - 1] + o2 / 2., ww_);

            k4 = DX_ * P_RHS(x_ + 2. * h, p[i - 1] + k3, dp[i - 1] + q3,
                             psi[i - 1] + s3, dpsi[i - 1] + r3,
                             omega[i - 1] + o3, ww_);
            q4 = DX_ * DP_RHS(x_ + 2. * h, p[i - 1] + k3, dp[i - 1] + q3,
                              psi[i - 1] + s3, dpsi[i - 1] + r3,
                              omega[i - 1] + o3, ww_);
            o4 = DX_ * OMEGA_RHS(x_ + 2. * h, p[i - 1] + k3, dp[i - 1] + q3,
                                 psi[i - 1] + s3, dpsi[i - 1] + r3,
                                 omega[i - 1] + o3, ww_);
            s4 = DX_ * PSI_RHS(x_ + 2. * h, p[i - 1] + k3, dp[i - 1] + q3,
                               psi[i - 1] + s3, dpsi[i - 1] + r3,
                               omega[i - 1] + o3, ww_);
            r4 = DX_ * DPSI_RHS(x_ + 2. * h, p[i - 1] + k3, dp[i - 1] + q3,
                                psi[i - 1] + s3, dpsi[i - 1] + r3,
                                omega[i - 1] + o3, ww_);

            index = i - 1;
            if (j == jmax)
            {
                index = i;
            }

            p[index] = p[i - 1] + (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
            dp[index] = dp[i - 1] + (q1 + 2. * q2 + 2. * q3 + q4) / 6.;
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
    double k1=0., k2=0., k3=0., k4=0., q1=0., q2=0., q3=0., q4=0., x_ = iter * dx, h,
                                           delta = (double)gridsize;
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
        o1 = DX_ * OMEGA_RHS(x_, p[i - 1], dp[i - 1], psi[i - 1], dpsi[i - 1],
                             omega[i - 1], ww_);
        s1 = DX_ * PSI_RHS(x_, p[i - 1], dp[i - 1], psi[i - 1], dpsi[i - 1],
                           omega[i - 1], ww_);
        r1 = DX_ * DPSI_RHS(x_, p[i - 1], dp[i - 1], psi[i - 1], dpsi[i - 1],
                            omega[i - 1], ww_);

        // k2 = dx*small_P_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] +
        // s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_); q2 = dx*DP_RHS(x_ +
        // h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] +
        // r1/2.,omega[i-1] + o1/2.,ww_);
        o2 = DX_ * OMEGA_RHS(x_ + h, p[i - 1] + k1 / 2., dp[i - 1] + q1 / 2.,
                             psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                             omega[i - 1] + o1 / 2., ww_);
        s2 = DX_ * PSI_RHS(x_ + h, p[i - 1] + k1 / 2., dp[i - 1] + q1 / 2.,
                           psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                           omega[i - 1] + o1 / 2., ww_);
        r2 = DX_ * DPSI_RHS(x_ + h, p[i - 1] + k1 / 2., dp[i - 1] + q1 / 2.,
                            psi[i - 1] + s1 / 2., dpsi[i - 1] + r1 / 2.,
                            omega[i - 1] + o1 / 2., ww_);

        // k3 = dx*small_P_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] +
        // s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_); q3 = dx*DP_RHS(x_ +
        // h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] +
        // r2/2.,omega[i-1] + o2/2.,ww_);
        o3 = DX_ * OMEGA_RHS(x_ + h, p[i - 1] + k2 / 2., dp[i - 1] + q2 / 2.,
                             psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                             omega[i - 1] + o2 / 2., ww_);
        s3 = DX_ * PSI_RHS(x_ + h, p[i - 1] + k2 / 2., dp[i - 1] + q2 / 2.,
                           psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                           omega[i - 1] + o2 / 2., ww_);
        r3 = DX_ * DPSI_RHS(x_ + h, p[i - 1] + k2 / 2., dp[i - 1] + q2 / 2.,
                            psi[i - 1] + s2 / 2., dpsi[i - 1] + r2 / 2.,
                            omega[i - 1] + o2 / 2., ww_);

        // k4 = dx*small_P_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] +
        // s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_); q4 = dx*DP_RHS(x_
        // + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] +
        // r3,omega[i-1] + o3,ww_);
        o4 = DX_ * OMEGA_RHS(x_ + 2. * h, p[i - 1] + k3, dp[i - 1] + q3,
                             psi[i - 1] + s3, dpsi[i - 1] + r3,
                             omega[i - 1] + o3, ww_);
        s4 = DX_ * PSI_RHS(x_ + 2. * h, p[i - 1] + k3, dp[i - 1] + q3,
                           psi[i - 1] + s3, dpsi[i - 1] + r3, omega[i - 1] + o3,
                           ww_);
        r4 = DX_ * DPSI_RHS(x_ + 2. * h, p[i - 1] + k3, dp[i - 1] + q3,
                            psi[i - 1] + s3, dpsi[i - 1] + r3,
                            omega[i - 1] + o3, ww_);

        p[i] = 0.;  // p[i-1] + (k1 + 2.*k2 + 2.*k3 + k4)/6.;
        dp[i] = 0.; // dp[i-1] + (q1 + 2.*q2 + 2.*q3 + q4)/6.;
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
        std::cout << "\33[30;41m"
                  << " Asymptotic Radius Too Small"
                  << "\x1B[0m" << std::endl;
        std::cout << x_ << std::endl;
    }
}

// these functions return the right hand side of the ode's
// small_P_RHS is valid for large redius when the scalar field is small.

double BosonStarSolution::small_P_RHS(const double x, const double P,
                                      const double DP, const double PSI,
                                      const double DPSI, const double OM,
                                      const double ww_)
{
    double RHS = -P * PSI * sqrt(DV(P) - ww_ / (OM * OM));
    return RHS;
}

double BosonStarSolution::P_RHS(const double x, const double P, const double DP,
                                const double PSI, const double DPSI,
                                const double OM, const double ww_)
{
    double RHS = DP;
    return RHS;
}

double BosonStarSolution::DP_RHS(const double x, const double P,
                                 const double DP, const double PSI,
                                 const double DPSI, const double OM,
                                 const double ww_)
{
    double r = ((x == 0.) ? eps : x);
    double DOM = OMEGA_RHS(x, P, DP, PSI, DPSI, OM, ww_);
    return P * PSI * PSI * (DV(P) - ww_ / (OM * OM)) -
           DP * (DOM / OM + DPSI / PSI + 2. / r);
}

double BosonStarSolution::PSI_RHS(const double x, const double P,
                                  const double DP, const double PSI,
                                  const double DPSI, const double OM,
                                  const double ww_)
{
    double RHS = DPSI;
    return RHS;
}

double BosonStarSolution::DPSI_RHS(const double x, const double P,
                                   const double DP, const double PSI,
                                   const double DPSI, const double OM,
                                   const double ww_)
{
    double r = ((x == 0.) ? eps : x);
    return 0.5 * DPSI * DPSI / PSI - 2. * DPSI / r -
           2. * M_PI * G * PSI *
               (PSI * PSI * V(P) + DP * DP +
                ww_ * P * P * PSI * PSI / (OM * OM));
}

double BosonStarSolution::OMEGA_RHS(const double x, const double P,
                                    const double DP, const double PSI,
                                    const double DPSI, const double OM,
                                    const double ww_)
{
    double r = ((x == 0.) ? eps : x);
    return (OM / (x * DPSI + PSI)) *
           (2. * M_PI * G * x * PSI *
                (DP * DP - PSI * PSI * V(P) +
                 ww_ * P * P * PSI * PSI / (OM * OM)) -
            DPSI - 0.5 * x * DPSI * DPSI / PSI);
}

// V is klein gordon potential and DV is its gradient. Depends on #define
// star_type at top
double BosonStarSolution::V(const double P)
{
    if (!solitonic)
    {
        return MM * P * P + 0.5 * lambda * P * P * P * P;
    }
    else
    {
        return MM * P * P * pow((1. - 2. * pow(P / sigma, 2)), 2);
    }
}
double BosonStarSolution::DV(const double P)
{
    if (!solitonic)
    {
        return MM + lambda * P * P;
    }
    else
    {
        return MM - 8. * MM * pow(P / sigma, 2) + 12. * MM * pow(P / sigma, 4);
    }
}

// counts how many times the function crosses the axis
int BosonStarSolution::crossings()
{
    int number = 0;
    for (int i = 0; i < gridsize - 1; ++i)
    {
        if (p[i] * p[i + 1] < 0)
        {
            number += 1;
        }
    }
    return number;
}

// 4th order error (cubic interpolation) for field. shouts if asked to fetch a
// value outside the ode solution
double BosonStarSolution::get_p_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 =
        ((iter == 0)
             ? p[1]
             : p[iter - 1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = p[iter];
    f3 = p[iter + 1];
    f4 = p[iter + 2];

    if (iter > gridsize - 3)
    {
        std::cout << "Requested Value outside BS initial data domain!"
                  << std::endl;
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

double BosonStarSolution::get_dp_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? dp[1] : dp[iter - 1]); // conditionl/ternary imposing
                                               // zero gradeint at r=0
    f2 = dp[iter];
    f3 = dp[iter + 1];
    f4 = dp[iter + 2];

    if (iter > gridsize - 3)
    {
        std::cout << "Requested Value outside BS initial data domain!"
                  << std::endl;
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
        std::cout << "Requested Value outside BS initial data domain!"
                  << std::endl;
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
        std::cout << "Requested Value outside BS initial data domain!"
                  << std::endl;
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
        std::cout << "Requested Value outside BS initial data domain!"
                  << std::endl;
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
        std::cout << "Requested Value outside BS initial data domain!"
                  << std::endl;
    }

    // do the cubic spline (for gradient now), from mathematica script written
    // by Robin (rc634@cam.ac.uk)
    interpolated_value =
        (1. / (24. * dx)) *
        ((f1 - 27. * f2 + 27. * f3 - f4) + 12. * a * (f1 - f2 - f3 + f4) -
         12. * a * a * (f1 - 3. * f2 + 3. * f3 - f4));
    return interpolated_value;
}

// returns the aspect mass, i.e. comparing the large radius metric to
// Schwarzschild and using a differential realtion to return M (NOT an ADM mass
// calc)
double BosonStarSolution::get_mass() const { return aspect_mass; }

// returns the eigenvalue
double BosonStarSolution::get_w() const { return sqrt(ww); }

double BosonStarSolution::get_r(const double frac) const
{
    if ((frac - 0.5) * (frac - 0.5) >= 0.25)
    {
        return -1.;
    }

    for (int i = 0; i < gridsize; ++i)
    {
        if (p[i] / p[0] < frac)
        {
            return radius_array[i];
        }
    }
    return -1.;
}

void BosonStarSolution::set_initialcondition_params(
    BosonStar_params_t m_params_BosonStar,
    Potential::params_t m_params_potential, const double max_r)
{
    gridsize = m_params_BosonStar.gridpoints;
    adaptive_buffer =
        0.; // gridsize/10; // numer of gridpoints to intergate more carefully
    p.resize(gridsize);            // scalar field modulus
    dp.resize(gridsize);           // scalar field modulus gradient
    psi.resize(gridsize);          // conformal factor
    dpsi.resize(gridsize);         // conformal factor gradient
    omega.resize(gridsize);        // lapse
    radius_array.resize(gridsize); // radius

    G = m_params_BosonStar.Newtons_constant;
    PC = m_params_BosonStar.central_amplitude_CSF;
    EIGEN = m_params_BosonStar.eigen;
    MM = m_params_potential.scalar_mass * m_params_potential.scalar_mass;
    lambda = m_params_potential.phi4_coeff;
    solitonic = m_params_potential.solitonic;
    sigma = m_params_potential.sigma_soliton;
    L = max_r * 1.05; // just to make sure the function domain is slightly
                      // larger than the required cube
    dx = L / (gridsize - 1);
}

void BosonStarSolution::shout() const
{
    std::cout << "Haliboombah!" << std::endl;
}

#endif /* BOSONSTARSOLUTION_IMPL_HPP_ */