/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARHELPERFUNCTION_HPP_
#define BOSONSTARHELPERFUNCTION_HPP_

#include "BosonStarParams.hpp"
#include "BosonStarSolver.hpp"
#include "ComplexPotential.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"

class BosonStarHelperFunction
{
  public:
    struct BS_3d_vars
    {
        double lapse;  // lapse
        double shiftx; // x-component of shift (all others are zero for a
                       // spherically symmetric BS)
        double phi_Re; // real part of \varphi
        double phi_Im; // imaginary part of \varphi
        double Pi_Re;  // real part of \Pi
        double Pi_Im;  // imaginary part of \Pi
        Tensor<2, double> KLL; // lower indices of extrinsic curvature
        Tensor<2, double> gLL; // lower indices of spatial ADM metric
    };

    // Write the solution variables to the vars object
    BS_3d_vars compute_star_vars(const Coordinates<double> &coords,
                                 double rapidity, double separation_offset,
                                 double impact_offset,
                                 const BosonStarSolver &sol,
                                 double phase_offset, bool is_antiboson) const
    {
        // Initialise
        BS_3d_vars vars;
        vars.lapse = 0.0;
        vars.shiftx = 0.0;
        vars.phi_Re = 0.0;
        vars.phi_Im = 0.0;
        vars.Pi_Re = 0.0;
        vars.Pi_Im = 0.0;
        FOR(i, j)
        {
            vars.KLL[i][j] = 0.0;
            vars.gLL[i][j] = 0.0;
        }

        // Positioning
        double c_ = cosh(rapidity);
        double s_ = sinh(rapidity);
        double v_ = tanh(rapidity);
        double t = (coords.x + separation_offset) * s_;
        double x = (coords.x + separation_offset) * c_;
        double y = coords.y + impact_offset;
        double z = coords.z;
        double r = sqrt(x * x + y * y + z * z);

        // Extract from 1d solver and interpolate onto the given coord
        double p_ = sol.get_A_interp(r);   // scalar amplitude
        double dp_ = sol.get_dA_interp(r); // derivative of the scalar amplitude
        double omega_ = sol.get_lapse_interp(r); // lapse
        double omega_prime_ =
            sol.get_dlapse_interp(r);        // derivative of the lapse
        double psi_ = sol.get_psi_interp(r); // conformal factor
        double psi_prime_ =
            sol.get_dpsi_interp(r); // derivative of the conformal factor
        double pc_os = psi_ * psi_ * c_ * c_ - omega_ * omega_ * s_ * s_;

        // BS frequency
        double w_ =
            is_antiboson ? -sol.get_BSfrequency() : sol.get_BSfrequency();

        // Gauge vars
        double lapse = omega_ * psi_ / sqrt(pc_os);
        vars.lapse = lapse;
        double phase_ = phase_offset * M_PI + w_ * t;
        double beta_x = s_ * c_ * (psi_ * psi_ - omega_ * omega_) / pc_os;
        vars.shiftx += beta_x;

        // Scalar vars
        vars.phi_Re += p_ * cos(phase_);
        vars.phi_Im += p_ * sin(phase_);
        vars.Pi_Re +=
            -(1. / lapse) * ((x / r) * (s_ - beta_x * c_) * dp_ * cos(phase_) -
                             w_ * (c_ - beta_x * s_) * p_ * sin(phase_));
        vars.Pi_Im +=
            -(1. / lapse) * ((x / r) * (s_ - beta_x * c_) * dp_ * sin(phase_) +
                             w_ * (c_ - beta_x * s_) * p_ * cos(phase_));

        // Metric components
        vars.gLL[0][0] = pc_os;
        vars.gLL[1][1] = psi_ * psi_;
        vars.gLL[2][2] = psi_ * psi_;

        // Extrinsic curvature components
        vars.KLL[2][2] = -lapse * s_ * x * psi_prime_ / (r * psi_);
        vars.KLL[1][1] = vars.KLL[2][2];
        vars.KLL[0][1] = lapse * c_ * s_ * (y / r) *
                         (psi_prime_ / psi_ - omega_prime_ / omega_);
        vars.KLL[0][2] = lapse * c_ * s_ * (z / r) *
                         (psi_prime_ / psi_ - omega_prime_ / omega_);
        vars.KLL[1][0] = vars.KLL[0][1];
        vars.KLL[2][0] = vars.KLL[0][2];
        vars.KLL[2][1] = 0.;
        vars.KLL[1][2] = 0.;
        vars.KLL[0][0] = lapse * (x / r) * s_ * c_ * c_ *
                         (psi_prime_ / psi_ - 2. * omega_prime_ / omega_ +
                          v_ * v_ * omega_ * omega_prime_ * pow(psi_, -2));

        return vars;
    }
};

#endif /* BOSONSTARHELPERFUNCTION_HPP_ */